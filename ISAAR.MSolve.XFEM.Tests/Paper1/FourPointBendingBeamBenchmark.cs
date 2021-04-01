using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.Paper1
{
    public class FourPointBendingBeamBenchmark
    {
        #region constants
        // Material and dimensions
        private const double fractureToughness = double.MaxValue;
        private const double h = 0.100, L = 0.440, thickness = 0.1;// m
        private const double v = 0.2, E = 3.5E7; // Pa

        // Boundary conditions
        private const double applicationArea = 0.020; // m
        private const double load = 20E3; // N
        private const double leftP = load / 11.0, offsetLeftP = 0.020;
        private const double rightP = load* 10.0 / 11.0, offsetRightP = 0.240;
        private const double offsetRoller = 0.200, offsetHinge = 0.420;

        // Crack
        private static readonly CartesianPoint crackMouth = new CartesianPoint(0.220, 0.100);
        private static readonly CartesianPoint crackTip = new CartesianPoint(0.220, 0.080);
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

        private readonly int numElementsX, numElementsY;
        private readonly int numSubdomainsX, numSubdomainsY;

        private readonly string plotDirectoryLsm, plotDirectorySubdomains;

        private readonly double tipEnrichmentRadius = 0.0;

        public FourPointBendingBeamBenchmark(int numElementsX, int numElementsY, int numSubdomainsX, int numSubdomainsY, 
            double growthLength, double tipEnrichmentRadius, double jIntegralRadiusOverElementSize, 
            int maxIterations, double heavisideTol, string plotDirectoryLsm, string plotDirectorySubdomains)
        {
            this.numElementsX = numElementsX;
            this.numElementsY = numElementsY;
            this.numSubdomainsX = numSubdomainsX;
            this.numSubdomainsY = numSubdomainsY;
            this.growthLength = growthLength;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.plotDirectoryLsm = plotDirectoryLsm;
            this.plotDirectorySubdomains = plotDirectorySubdomains;
            this.MaxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description. Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public TrackingExteriorCrackLsm Crack { get; private set; }

        public double FractureToughness => double.MaxValue;

        public int MaxIterations { get; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public XModel Model { get; private set; }

        public string Name => "4-point bending beam benchmark";

        public void InitializeModel()
        {
            Model = new XModel();
            CreateModel();
            InitializeCrack();
        }

        public void CreateModel()
        {
            var builder = new Uniform2DXModelBuilder();
            builder.DomainLengthX = L;
            builder.DomainLengthY = h;
            builder.NumSubdomainsX = numSubdomainsX;
            builder.NumSubdomainsY = numSubdomainsY;
            builder.NumTotalElementsX = numElementsX;
            builder.NumTotalElementsY = numElementsY;
            builder.PlaneStress = true;
            builder.YoungModulus = E;
            builder.PoissonRatio = v;
            builder.Thickness = thickness;

            Model = builder.BuildModel();

            // Boundary conditions
            double meshTol = 1E-8;
            double start, end;
            double load;
            IEnumerable<XNode> nodes;
            IEnumerable<XNode> bottomNodes = Model.Nodes.Values.Where(n => Math.Abs(n.Y) <= meshTol);
            IEnumerable<XNode> topNodes = Model.Nodes.Values.Where(n => Math.Abs(n.Y - h) <= meshTol);

            // Left load
            start = offsetLeftP - applicationArea / 2.0 - meshTol;
            end = offsetLeftP + applicationArea / 2.0 + meshTol;
            nodes = bottomNodes.Where(n => (n.X >= start) && (n.X <= end));
            load = leftP / nodes.Count();
            foreach (XNode node in nodes) Model.NodalLoads.Add(new NodalLoad(node, StructuralDof.TranslationY, load));

            // Right load
            start = offsetRightP - applicationArea / 2.0 - meshTol;
            end = offsetRightP + applicationArea / 2.0 + meshTol;
            nodes = bottomNodes.Where(n => (n.X >= start) && (n.X <= end));
            load = rightP / nodes.Count();
            foreach (XNode node in nodes) Model.NodalLoads.Add(new NodalLoad(node, StructuralDof.TranslationY, load));

            // Roller BC
            start = offsetRoller - applicationArea / 2.0 - meshTol;
            end = offsetRoller + applicationArea / 2.0 + meshTol;
            nodes = topNodes.Where(n => (n.X >= start) && (n.X <= end));
            foreach (XNode node in nodes) node.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY, Amount = 0.0 });

            // Hinge BC
            start = offsetHinge - applicationArea / 2.0 - meshTol;
            end = offsetHinge + applicationArea / 2.0 + meshTol;
            nodes = topNodes.Where(n => (n.X >= start) && (n.X <= end));
            foreach (XNode node in nodes)
            {
                node.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                node.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            }
        }

        public void InitializeCrack()
        {
            var globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(0, E, v);
            IPropagator propagator = new Propagator(Model.Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            var initialCrack = new PolyLine2D(crackMouth, crackTip);

            var lsmCrack = new TrackingExteriorCrackLsm(propagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol),
                new SignFunction2D());
            lsmCrack.Mesh = Model.Mesh;

            // Logging         
            if (plotDirectoryLsm != null)
            {
                lsmCrack.EnrichmentLogger = new EnrichmentLogger(Model, lsmCrack, plotDirectoryLsm);
                lsmCrack.LevelSetLogger = new LevelSetLogger(Model, lsmCrack, plotDirectoryLsm);
                lsmCrack.LevelSetComparer = new PreviousLevelSetComparer(lsmCrack, plotDirectoryLsm);
            }

            // Mesh geometry interaction
            lsmCrack.InitializeGeometry(initialCrack);
            this.Crack = lsmCrack;
        }

        public class Builder //: IBenchmarkBuilder
        {
            private readonly int numElementsX, numElementsY, numSubdomainsX, numSubdomainsY;

            public Builder(int numElementsX, int numElementsY, int numSubdomainsX, int numSubdomainsY)
            {
                this.numElementsX = numElementsX;
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
            public double GrowthLength { get; set; } = 0.01; //m

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string LsmPlotDirectory { get; set; } = null;

            public string SubdomainPlotDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = 8;

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public FourPointBendingBeamBenchmark BuildBenchmark()
            {
                return new FourPointBendingBeamBenchmark(numElementsX,numElementsY, numSubdomainsX, numSubdomainsY, GrowthLength,
                    TipEnrichmentRadius, JintegralRadiusOverElementSize, MaxIterations, HeavisideEnrichmentTolerance,
                    LsmPlotDirectory, SubdomainPlotDirectory);
            }
        }
    }
}
