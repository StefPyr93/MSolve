using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.GMSH;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests
{
    public class DcbBenchmarkBelytschkoMpi //: IBenchmark
    {
        #region constants
        ///// <summary>
        ///// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        ///// </summary>
        //private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
        //    HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        public static readonly double h = 3.94, L = 3 * h;// in
        private static readonly double v = 0.3, E = 3e7; // psi=lbs/in^2
        private static readonly double load = 197; // lbs
        private static readonly double a = 3.95, da = 0.5; // in 
        private static readonly double dTheta = 5.71 * Math.PI / 180; // initial crack angle
        #endregion

        private readonly double heavisideTol;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        private readonly int numElementsY;
        private readonly int numSubdomainsX;
        private readonly int numSubdomainsY;

        private readonly double tipEnrichmentRadius = 0.0;

        private readonly ProcessDistribution procs;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        public DcbBenchmarkBelytschkoMpi(ProcessDistribution procs, int numElementsY, int numSubdomainsX, int numSubdomainsY,
            double growthLength, double tipEnrichmentRadius, double jIntegralRadiusOverElementSize, int maxIterations, 
            double heavisideTol)
        {
            this.procs = procs;
            this.numElementsY = numElementsY;
            this.numSubdomainsX = numSubdomainsX;
            this.numSubdomainsY = numSubdomainsY;
            this.growthLength = growthLength;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.MaxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        public static IXFiniteElementFactory ElementFactory
        {
            get
            {
                var modelBuilder = new Uniform2DXModelBuilder();
                var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(0, E, v);
                var elementFactory = new XContinuumElement2DFactory(modelBuilder.QuadratureForStiffness,
                    modelBuilder.QuadratureForJintegral, material);
                return elementFactory;
            }
        }

        /// <summary>
        /// The crack geometry description. Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public ICrackDescriptionMpi Crack { get; private set; }

        public double FractureToughness => fractureToughness;

        //public IReadOnlyList<double> GrowthAngles { get; private set; }
        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        public int MaxIterations { get; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public IXModelMpi Model { get; private set; }

        public string Name { get { return "Fillet benchmark"; } }

        //public string PlotDirectory { get { return lsmPlotDirectory; } }

        public TipAdaptivePartitioner Partitioner { get; set; } // Refactor its injection

        public void InitializeModel()
        {
            CreateModel();
            InitializeCrack();
        }

        public void CreateModel()
        {
            Func<XModel> createModel = () =>
            {
                var builder = new Uniform2DXModelBuilder();
                builder.DomainLengthX = L;
                builder.DomainLengthY = h;
                builder.NumSubdomainsX = numSubdomainsX;
                builder.NumSubdomainsY = numSubdomainsY;
                builder.NumTotalElementsX = 3 * numElementsY;
                builder.NumTotalElementsY = numElementsY;
                builder.PlaneStress = false;
                builder.YoungModulus = E;
                builder.PoissonRatio = v;
                builder.PrescribeDisplacement(Uniform2DXModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationX, 0.0);
                builder.PrescribeDisplacement(Uniform2DXModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 0.0);
                builder.DistributeLoadAtNodes(Uniform2DXModelBuilder.BoundaryRegion.UpperLeftCorner, StructuralDof.TranslationY, load);
                builder.DistributeLoadAtNodes(Uniform2DXModelBuilder.BoundaryRegion.LowerLeftCorner, StructuralDof.TranslationY, -load);

                return builder.BuildModel();
            };

            //Model = new XModelMpiCentralized(procs, createModel, ElementFactory);
            Model = new XModelMpiRedundant(procs, createModel);
        }

        public void InitializeCrack()
        {
            Func<TrackingExteriorCrackLsm> createCrack = () =>
            {
                var model = Model.RawModel;
                var mesh = Model.RawModel.Mesh;
                var globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(0, E, v);
                IPropagator propagator = new Propagator(Model.RawModel.Mesh, jIntegralRadiusOverElementSize,
                    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

                var crackMouth = new CartesianPoint(0.0, h / 2);
                var crackKink = new CartesianPoint(a, h / 2);
                var initialCrack = new PolyLine2D(crackMouth, crackKink);
                initialCrack.UpdateGeometry(-dTheta, da);
                //var crackTip = new CartesianPoint(a + da * Math.Cos(dTheta), h/2 - da * Math.Sin(dTheta));

                TrackingExteriorCrackLsm lsmCrack = new TrackingExteriorCrackLsm(propagator, tipEnrichmentRadius, 
                    new RelativeAreaResolver(heavisideTol), new SignFunction2D());
                lsmCrack.Mesh = Model.RawModel.Mesh;

                // Mesh geometry interaction
                lsmCrack.InitializeGeometry(initialCrack);
                //lsmCrack.UpdateGeometry(-dTheta, da);

                return lsmCrack;
            };

            //this.Crack = new TrackingExteriorCrackLsmMpiCentralized(procs, createCrack);
            this.Crack = new TrackingExteriorCrackLsmMpiRedundant(procs, createCrack);
            Model.DofSerializer = new EnrichedDofSerializer(this.Crack);
        }

        public class Builder //: IBenchmarkBuilder
        {
            private readonly int numElementsY;
            private readonly int numSubdomainsX;
            private readonly int numSubdomainsY;
            private readonly ProcessDistribution procs;

            public Builder(ProcessDistribution procs, int numElementsY, int numSubdomainsX, int numSubdomainsY)
            {
                this.procs = procs;
                this.numElementsY = numElementsY;
                this.numSubdomainsX = numSubdomainsX;
                this.numSubdomainsY = numSubdomainsY;
            }

            /// <summary>
            /// A node that lies in the positive halfplane defined by the body level set of a crack, will be enriched with 
            /// heaviside enrichment if Apos/(Apos+Aneg) &gt; <see cref="HeavisideEnrichmentTolerance"/> where Apos, Aneg 
            /// are the subareas of its nodal support  in the positive and negative halfplanes respectively. Similarly a
            /// node in the negative halfplane will be enriched if Aneg/(Apos+Aneg) &gt; 
            /// <see cref="HeavisideEnrichmentTolerance"/>.
            /// </summary>
            public double HeavisideEnrichmentTolerance { get; set; } = 0.0001;

            /// <summary>
            /// Must be sufficiently larger than the element size.
            /// </summary>
            public double GrowthLength { get; set; } = 0.3; //in

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = 8; //TODO: After that I noticed very weird behaviour

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public DcbBenchmarkBelytschkoMpi BuildBenchmark()
            {
                return new DcbBenchmarkBelytschkoMpi(procs, numElementsY, numSubdomainsX, numSubdomainsY, GrowthLength, 
                    TipEnrichmentRadius, JintegralRadiusOverElementSize, MaxIterations, HeavisideEnrichmentTolerance);
            }
        }
    }
}
