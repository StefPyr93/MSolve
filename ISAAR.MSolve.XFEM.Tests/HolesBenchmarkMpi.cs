using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;
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
    public class HolesBenchmarkMpi
    {
        public enum BoundaryConditions
        {
            BottomConstrainXY_TopDisplacementY, BottomConstrainXY_TopConstrainXDisplacementY,
            BottomConstrainY_TopDisplacementY, BottomDisplacementY_TopDisplacementY,
            BottomConstrainXDisplacementY_TopConstrainXDisplacementY
        }

        #region constants
        /// <summary>
        /// Thickness of the whole domain
        /// </summary>
        private const double t = 1.0; // mm

        /// <summary>
        /// Young's modulus
        /// </summary>
        private const double E = 2e5; // N/mm^2

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// Fracture toughness
        /// </summary>
        private const double KIc = double.MaxValue; // N.mm^(3/2) actually 1500

        /// <summary>
        /// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        /// </summary>
        private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
            HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(0, E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        /// <summary>
        /// Prescribed displacement
        /// </summary>
        private const double displacement = 0.1; //mm

        // Geometry
        private const double minX = 0.0, minY = 0.0, maxX = 20.0, maxY = 10.0; //mm
        private const double holeRadius = 2.0, leftHoleX = 3.0, leftHoleY = 7.0, rightHoleX = 17.0, rightHoleY = 3.0; //mm
        private const double initialCrackLength = 1.0;
        private const double leftCrackMouthX = 0.0, leftCrackMouthY = 2.85, leftCrackTipX = minX + initialCrackLength, leftCrackTipY = 2.85; //mm
        private const double rightCrackMouthX = 20.0, rightCrackMouthY = 7.15, rightCrackTipX = maxX - initialCrackLength, rightCrackTipY = 7.15; //mm
        
        private const int subdomainID = 0;
        #endregion

        private readonly BoundaryConditions bc;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;
        private readonly double heavisideTol;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;
        private readonly string meshPath;
        private readonly ProcessDistribution procs;
        private readonly double tipEnrichmentRadius;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private HolesBenchmarkMpi(ProcessDistribution procs, string meshPath, double growthLength, BoundaryConditions bc, 
            double jIntegralRadiusOverElementSize, double tipEnrichmentRadius, int maxIterations, double heavisideTol)
        {
            this.procs = procs;
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.bc = bc;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.MaxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescriptionMpi Crack { get; private set; }

        public double FractureToughness => fractureToughness;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        public int MaxIterations { get; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public IXModelMpi Model { get; private set; }

        public string Name { get { return "Twin Holes benchmark"; } }

        public ISingleCrack LeftCrack => Crack.SingleCracks[0];
        public ISingleCrack RightCrack => Crack.SingleCracks[1];

        public TipAdaptivePartitioner Partitioner { get; set; } // Refactor its injection

        public void InitializeModel(int numSubdomains)
        {
            CreateModel(numSubdomains);
            InitializeCrack();
        }

        public void WritePropagation()
        {
            // Write crack path
            Console.WriteLine("Left crack path:");
            foreach (var point in LeftCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
            Console.WriteLine("Right crack path:");
            foreach (var point in RightCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
        }

        private void ApplyBoundaryConditions(XModel model)
        {
            double meshTol = 1E-6;
            XNode leftTopCorner = model.Nodes.Values.Where(
                node => (Math.Abs(node.X - minX) <= meshTol) && (Math.Abs(node.Y - maxY) <= meshTol)).First();
            XNode rightBottomCorner = model.Nodes.Values.Where(
                node => (Math.Abs(node.X - maxX) <= meshTol) && (Math.Abs(node.Y - minY) <= meshTol)).First();
            XNode[] bottomNodes = model.Nodes.Values.Where(node => Math.Abs(node.Y - minY) <= meshTol).ToArray();
            XNode[] topNodes = model.Nodes.Values.Where(node => Math.Abs(node.Y - maxY) <= meshTol).ToArray();

            if (bc == BoundaryConditions.BottomConstrainXY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainXY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
                leftTopCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                rightBottomCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            }
            else if (bc == BoundaryConditions.BottomDisplacementY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.5 * displacement });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.5 * displacement });
                }
                leftTopCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                rightBottomCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            }
            else if (bc == BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.5 * displacement });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.5 * displacement });
                }
            }
            else throw new Exception("Shouldn't have been reached.");
        }

        private void CreateModel(int numSubdomains)
        {
            // Element type selection
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(0, E, v);
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, 2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, 2));
            var jIntegration =
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, 4);
            var elementFactory = new XContinuumElement2DFactory(integration, jIntegration, material);

            Func<XModel> createModel = () =>
            {
                var model = new XModel();
                model.Subdomains.Add(subdomainID, new XSubdomain(subdomainID));

                // Mesh generation
                var reader = new GmshReader<XNode>(meshPath);
                (IReadOnlyList<XNode> nodes, IReadOnlyList<CellConnectivity<XNode>> elementConnectivities) = reader.CreateMesh(
                    (id, x, y, z) => new XNode(id, x, y, z));

                // Nodes
                foreach (XNode node in nodes) model.Nodes.Add(node.ID, node);

                // Elements
                var cells = new XContinuumElement2D[elementConnectivities.Count];
                for (int e = 0; e < cells.Length; ++e)
                {
                    XContinuumElement2D element =
                        elementFactory.CreateElement(e, CellType.Quad4, elementConnectivities[e].Vertices);
                    cells[e] = element;
                    model.Elements.Add(e, element);
                    model.Subdomains[subdomainID].Elements.Add(e, model.Elements[e]);
                }

                // Mesh usable for crack-mesh interaction
                var boundary = new HolesBoundary();
                model.Boundary = boundary;
                model.Mesh =
                    new BidirectionalMesh2D<XNode, XContinuumElement2D>(model.Nodes.Values.ToArray(), cells, boundary);

                // Apply boundary conditions
                ApplyBoundaryConditions(model);

                // Partition into subdomain
                PartitionMesh(numSubdomains, model, model.Mesh);

                return model;
            };

            //Model = new XModelMpiCentralized(procs, createModel, elementFactory);
            Model = new XModelMpiRedundant(procs, createModel);
        }

        private void InitializeCrack()
        {
            var createSingleCracks = new Func<TrackingExteriorCrackLsm>[2];
            createSingleCracks[0] = () =>
            {
                // Left crack
                IPropagator leftPropagator = new Propagator(Model.RawModel.Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

                var initialLeftCrack = new PolyLine2D(new CartesianPoint(leftCrackMouthX, leftCrackMouthY),
                    new CartesianPoint(leftCrackTipX, leftCrackTipY));
                TrackingExteriorCrackLsm leftCrack = new TrackingExteriorCrackLsm(leftPropagator, tipEnrichmentRadius,
                    new RelativeAreaResolver(heavisideTol), new SignFunction2D());
                leftCrack.Mesh = Model.RawModel.Mesh;
                leftCrack.InitializeGeometry(initialLeftCrack);

                return leftCrack;
            };

            createSingleCracks[1] = () =>
            { 
                // Right crack
                IPropagator rightPropagator = new Propagator(Model.RawModel.Mesh, jIntegralRadiusOverElementSize,
                    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

                var initialRightCrack = new PolyLine2D(new CartesianPoint(rightCrackMouthX, rightCrackMouthY),
                    new CartesianPoint(rightCrackTipX, rightCrackTipY));
                TrackingExteriorCrackLsm rightCrack = new TrackingExteriorCrackLsm(rightPropagator, tipEnrichmentRadius,
                    new RelativeAreaResolver(heavisideTol), new SignFunction2D());
                rightCrack.Mesh = Model.RawModel.Mesh;
                rightCrack.InitializeGeometry(initialRightCrack);

                return rightCrack;
            };

            var singleCracks = new TrackingExteriorCrackLsm[2];
            singleCracks[0] = createSingleCracks[0]();
            singleCracks[1] = createSingleCracks[1]();

            //this.Crack = new MultipleCracksDisjointMpiCentralized(procs, singleCracks);
            this.Crack = new MultipleCracksDisjointMpiRedundant(procs, singleCracks);
            Model.DofSerializer = new EnrichedDofSerializer(this.Crack);
        }

        private static void PartitionMesh(int numSubdomains, XModel model, IMesh2D<XNode, XContinuumElement2D> mesh)
        {
            // Define subdomain boundaries
            double tol = 1E-13;
            var regions = new Dictionary<int, IRegion2D>();
            if (numSubdomains == 10)
            {
                double xMin = 0.0, yMin = 0.0, xMax = 20.0, yMax = 10.0;
                double x1 = 2.75, x2 = 3.5, x3 = 5.0, x4 = 8.0, x5 = 12.2, x6 = 15.5, x7 = 16.5, x8 = 16.75;
                double y1 = 1.85, y2 = 5.0, y3 = 8.85;

                // Left crack regions: 
                var region0 = new RectangularRegion2D(xMin, yMin, x2, y2, tol);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[0] = region0;

                var region1 = new RectangularRegion2D(x2, yMin, x4, y2, tol);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = region1;

                var region2 = new RectangularRegion2D(x4, yMin, x5, y2, tol);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2] = region2;

                var region3Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x5, yMin), new CartesianPoint(x6, yMin), new CartesianPoint(x6, y1),
                    new CartesianPoint(x8, y2), new CartesianPoint(x5, y2)
                };
                var region3Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region3Vertices[1], region3Vertices[2]),
                    new LineSegment2D(region3Vertices[2], region3Vertices[3]),
                    new LineSegment2D(region3Vertices[3], region3Vertices[4]),
                    new LineSegment2D(region3Vertices[4], region3Vertices[1])
                };
                regions[3] = new PolygonalRegion2D(region3Vertices, region3Boundaries);

                var region4Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x6, yMin), new CartesianPoint(xMax, yMin), new CartesianPoint(xMax, y2),
                    new CartesianPoint(x8, y2), new CartesianPoint(x6, y1)
                };
                var region4Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region4Vertices[2], region4Vertices[3]),
                    new LineSegment2D(region4Vertices[3], region4Vertices[4]),
                    new LineSegment2D(region4Vertices[4], region4Vertices[0])
                };
                regions[4] = new PolygonalRegion2D(region4Vertices, region4Boundaries);

                // Right crack regions:
                var region5Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(xMin, y2), new CartesianPoint(x1, y2), new CartesianPoint(x3, y3),
                    new CartesianPoint(x3, yMax), new CartesianPoint(xMin, yMax)
                };
                var region5Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region5Vertices[0], region5Vertices[1]),
                    new LineSegment2D(region5Vertices[1], region5Vertices[2]),
                    new LineSegment2D(region5Vertices[2], region5Vertices[3]),
                };
                regions[5] = new PolygonalRegion2D(region5Vertices, region5Boundaries);

                var region6Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x1, y2), new CartesianPoint(x4, y2), new CartesianPoint(x4, yMax),
                    new CartesianPoint(x3, yMax), new CartesianPoint(x3, y3)
                };
                var region6Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region6Vertices[0], region6Vertices[1]),
                    new LineSegment2D(region6Vertices[1], region6Vertices[2]),
                    new LineSegment2D(region6Vertices[3], region6Vertices[4]),
                    new LineSegment2D(region6Vertices[4], region6Vertices[0]),
                };
                regions[6] = new PolygonalRegion2D(region6Vertices, region6Boundaries);

                var region7 = new RectangularRegion2D(x4, y2, x5, yMax, tol);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[7] = region7;

                var region8 = new RectangularRegion2D(x5, y2, x7, yMax, tol);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[8] = region8;

                var region9 = new RectangularRegion2D(x7, y2, xMax, yMax, tol);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[9] = region9;
            }
            else if (numSubdomains == 15)
            {
                double xMin = 0.0, yMin = 0.0, xMax = 20.0, yMax = 10.0;
                double x1 = 2.75, x2 = 3.5, x3 = 5.0, x4 = 8.0, x5 = 12.2, x6 = 15.5, x7 = 16.5, x8 = 16.75;
                double y1 = 1.85, y2 = 3.25, y3 = 6.5, y4 = 8.85;

                // Bottom row: 
                var region0 = new RectangularRegion2D(xMin, yMin, x2, y2, tol);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[0] = region0;

                var region1 = new RectangularRegion2D(x2, yMin, x4, y2, tol);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = region1;

                var region2 = new RectangularRegion2D(x4, yMin, x5, y2, tol);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2] = region2;

                var region3Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x5, yMin), new CartesianPoint(x6, yMin), new CartesianPoint(x6, y1),
                    new CartesianPoint(x8, y2), new CartesianPoint(x5, y2)
                };
                var region3Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region3Vertices[1], region3Vertices[2]),
                    new LineSegment2D(region3Vertices[2], region3Vertices[3]),
                    new LineSegment2D(region3Vertices[3], region3Vertices[4]),
                    new LineSegment2D(region3Vertices[4], region3Vertices[1])
                };
                regions[3] = new PolygonalRegion2D(region3Vertices, region3Boundaries);

                var region4Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x6, yMin), new CartesianPoint(xMax, yMin), new CartesianPoint(xMax, y2),
                    new CartesianPoint(x8, y2), new CartesianPoint(x6, y1)
                };
                var region4Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region4Vertices[2], region4Vertices[3]),
                    new LineSegment2D(region4Vertices[3], region4Vertices[4]),
                    new LineSegment2D(region4Vertices[4], region4Vertices[0])
                };
                regions[4] = new PolygonalRegion2D(region4Vertices, region4Boundaries);

                // Middle row:
                var region5 = new RectangularRegion2D(xMin, y2, x2, y3, tol);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[5] = region5;

                var region6 = new RectangularRegion2D(x2, y2, x4, y3, tol);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[6] = region6;

                var region7 = new RectangularRegion2D(x4, y2, x5, y3, tol);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[7] = region7;

                var region8 = new RectangularRegion2D(x5, y2, x7, y3, tol);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[8] = region8;

                var region9 = new RectangularRegion2D(x7, y2, xMax, y3, tol);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[9] = region9;

                // Top row:
                var region10Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(xMin, y3), new CartesianPoint(x1, y3), new CartesianPoint(x3, y4),
                    new CartesianPoint(x3, yMax), new CartesianPoint(xMin, yMax)
                };
                var region10Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region10Vertices[0], region10Vertices[1]),
                    new LineSegment2D(region10Vertices[1], region10Vertices[2]),
                    new LineSegment2D(region10Vertices[2], region10Vertices[3]),
                };
                regions[10] = new PolygonalRegion2D(region10Vertices, region10Boundaries);

                var region11Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x1, y3), new CartesianPoint(x4, y3), new CartesianPoint(x4, yMax),
                    new CartesianPoint(x3, yMax), new CartesianPoint(x3, y4)
                };
                var region11Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region11Vertices[0], region11Vertices[1]),
                    new LineSegment2D(region11Vertices[1], region11Vertices[2]),
                    new LineSegment2D(region11Vertices[3], region11Vertices[4]),
                    new LineSegment2D(region11Vertices[4], region11Vertices[0]),
                };
                regions[11] = new PolygonalRegion2D(region11Vertices, region11Boundaries);

                var region12 = new RectangularRegion2D(x4, y3, x5, yMax, tol);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[12] = region12;

                var region13 = new RectangularRegion2D(x5, y3, x7, yMax, tol);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[13] = region13;

                var region14 = new RectangularRegion2D(x7, y3, xMax, yMax, tol);
                region14.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region14.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[14] = region14;
            }
            else
            {
                throw new NotImplementedException();
            }

            //Console.WriteLine($"Before partitioning: {benchmark.Model.NumSubdomains} subdomains");

            // Partition mesh into subdomains
            var regionsGeneral = new Dictionary<int, IRegion2D>();
            foreach (var pair in regions) regionsGeneral[pair.Key] = pair.Value;
            var partitioner = new GuidedPartioner2D<XNode, XContinuumElement2D>(mesh, regionsGeneral);
            Dictionary<int, List<XContinuumElement2D>> elementsOfSubdomains = partitioner.CreateSubdomains();

            // Replace the single subdomain that was already created with the ones from the mesh partition
            model.Subdomains.Clear();
            foreach (int subdomainID in elementsOfSubdomains.Keys)
            {
                model.Subdomains.Add(subdomainID, new XSubdomain(subdomainID));
                foreach (XContinuumElement2D element in elementsOfSubdomains[subdomainID])
                {
                    model.Subdomains[subdomainID].Elements[element.ID] = element;
                }
            }

            //Console.WriteLine($"After partitioning: {benchmark.Model.NumSubdomains} subdomains");
        }

        public class Builder //: IBenchmarkBuilder
        {
            private readonly double growthLength;
            private readonly string meshPath;
            private readonly ProcessDistribution procs;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="meshPath">The absolute path of the mesh file.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            /// <param name="timingDirectory">The absolute path of the file where slover timing will be written.</param>
            public Builder(ProcessDistribution procs, string meshPath, double growthLength)
            {
                this.procs = procs;
                this.growthLength = growthLength;
                this.meshPath = meshPath;
            }

            public BoundaryConditions BC { get; set; } = BoundaryConditions.BottomConstrainXY_TopDisplacementY;

            /// <summary>
            /// A node that lies in the positive halfplane defined by the body level set of a crack, will be enriched with 
            /// heaviside enrichment if Apos/(Apos+Aneg) &gt; <see cref="HeavisideEnrichmentTolerance"/> where Apos, Aneg 
            /// are the subareas of its nodal support  in the positive and negative halfplanes respectively. Similarly a
            /// node in the negative halfplane will be enriched if Aneg/(Apos+Aneg) &gt; 
            /// <see cref="HeavisideEnrichmentTolerance"/>.
            /// </summary>
            public double HeavisideEnrichmentTolerance { get; set; } = 0.0001;

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public HolesBenchmarkMpi BuildBenchmark()
            {
                return new HolesBenchmarkMpi(procs, meshPath, growthLength, BC, JintegralRadiusOverElementSize, 
                    TipEnrichmentRadius, MaxIterations, HeavisideEnrichmentTolerance);
            }
        }

        private class HolesBoundary : IDomain2DBoundary
        {
            public HolesBoundary()
            {
            }

            public bool IsInside(CartesianPoint point)
            {
                // Shapes
                var rectHull = new Rectangular2DBoundary(minX, maxX, minY, maxY);
                var leftCircle = new Circle2D(new CartesianPoint(leftHoleX, leftHoleY), holeRadius);
                var rightCircle = new Circle2D(new CartesianPoint(rightHoleX, rightHoleY), holeRadius);

                // Intrnal points lie inside the rectangle, but outside the circular holes.
                if (rectHull.IsInside(point))
                {
                    if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                    if (rightCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                }
                return false;
            }
        }
    }
}