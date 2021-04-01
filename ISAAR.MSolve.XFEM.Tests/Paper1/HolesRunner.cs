using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using static ISAAR.MSolve.XFEM.Tests.COMPDYN2019.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Paper1
{
    public class HolesRunner
    {
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_4442.msh";
        private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_8000.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_13738.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_22666.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_29052.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_60000.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_100000.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_189740.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Mesh\holes_357324.msh";
        private const string leftCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Plots\Left";
        private const string rightCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Plots\Right";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\Holes\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\solver_log.txt";

        public static void RunTest()
        {
            bool plotLSM = false;
            bool plotSubdomains = false;
            bool reanalysis = true;


            //HolesBenchmark benchmarkSub1;
            HolesBenchmark benchmarkSub10;
            //HolesBenchmark benchmarkSub15;

            // FETI-DP 10 subdomains
            benchmarkSub10 = CreateMultiSubdomainBenchmark(10, plotLSM);
            ISolverMpi solverFetiDP = DefineSolver(benchmarkSub10);
            RunCrackPropagationAnalysis(benchmarkSub10, solverFetiDP, plotSubdomains, reanalysis);
            //PlotSubdomains(benchmarkSub10, solverFetiDP);
            //Console.WriteLine("Uncracked analysis, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub10.Model, solverFetiDP));
            //Console.WriteLine("Cracked analysis only 1 step, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub10.Model, benchmarkSub10.Crack, solverFetiDP));
            //Console.WriteLine("FETI-DP, 10 subdomains: ");

            // FETI-DP 15 subdomains
            //benchmarkSub15 = CreateMultiSubdomainBenchmark(15);
            //ISolver solverFetiDP = DefineSolver(benchmarkSub15, SolverType.FetiDP);
            //RunCrackPropagationAnalysis(benchmarkSub15, solverFetiDP);

            Console.Write("\nEnd");
        }
        private static HolesBenchmark CreateSingleSubdomainBenchmark(bool plotLSM)
        {
            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.

            var builder = new HolesBenchmark.Builder(meshPath, growthLength);
            if (plotLSM)
            {
                builder.LeftLsmPlotDirectory = leftCrackPlotDirectory;
                builder.RightLsmPlotDirectory = rightCrackPlotDirectory;
            }

            builder.HeavisideEnrichmentTolerance = 0.12;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // If you modify the following two parameters significantly, then you will need to redefine which nodes are expected 
            // to be enriched.
            builder.TipEnrichmentRadius = 0.5;
            builder.BC = HolesBenchmark.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;

            builder.MaxIterations = 12;

            HolesBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static HolesBenchmark CreateMultiSubdomainBenchmark(int numSubdomains, bool plotLSM)
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

            // Create the model, crack & mesh with only 1 subdomain
            HolesBenchmark benchmark = CreateSingleSubdomainBenchmark(plotLSM);

            // Partition mesh into subdomains
            var regionsGeneral = new Dictionary<int, IRegion2D>();
            foreach (var pair in regions) regionsGeneral[pair.Key] = pair.Value;
            var partitioner = new GuidedPartioner2D<XNode, XContinuumElement2D>(benchmark.Mesh, regionsGeneral);
            Dictionary<int, List<XContinuumElement2D>> elementsOfSubdomains = partitioner.CreateSubdomains();

            // Replace the single subdomain that was already created with the ones from the mesh partition
            benchmark.Model.Subdomains.Clear();
            foreach (int subdomainID in elementsOfSubdomains.Keys)
            {
                benchmark.Model.Subdomains.Add(subdomainID, new XSubdomain(subdomainID));
                foreach (XContinuumElement2D element in elementsOfSubdomains[subdomainID])
                {
                    benchmark.Model.Subdomains[subdomainID].Elements[element.ID] = element;
                }
            }

            return benchmark;
        }

        private static ISolverMpi DefineSolver(HolesBenchmark benchmark)
        {
            benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
            //Dictionary<ISubdomain, HashSet<INode>> cornerNodes = null;
            Func<ISubdomain, HashSet<INode>> getCornerNodes = null;

            if (benchmark.Model.Subdomains.Count == 10 || benchmark.Model.Subdomains.Count == 15)
            {
                //cornerNodes = FindCornerNodesFromCrosspoints2D(benchmark.Model);
                getCornerNodes = (sub) => new HashSet<INode>(sub.EnumerateNodes().Where(node => node.Multiplicity > 2));
            }
            else
            {
                throw new NotImplementedException();
            }

            // Must also specify corner nodes
            var cornerNodeSelection = new CrackedFetiDPCornerNodesSerial(benchmark.Model, benchmark.Crack, getCornerNodes);
            var reordering = new OrderingAmdCSparseNet();
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            var builder = new FetiDPSolverSerial.Builder(fetiMatrices);
            //builder.Preconditioning = new LumpedPreconditioning();
            //builder.Preconditioning = new DiagonalDirichletPreconditioning();
            builder.Preconditioning = new DirichletPreconditioning();
            builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
            return builder.Build(benchmark.Model, cornerNodeSelection);
        }

        private static void RunCrackPropagationAnalysis(HolesBenchmark benchmark, ISolverMpi solver, 
            bool plotSubdomains, bool reanalysis)
        {
            var analyzer = new QuasiStaticCrackPropagationAnalyzerSerial(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, reanalysis, benchmark.Partitioner);

            // Subdomain plots
            if (plotSubdomains)
            {
                if (solver is FetiDPSolverSerial fetiDP)
                {
                    analyzer.DDLogger = new DomainDecompositionLoggerFetiDP(subdomainPlotDirectory, fetiDP.CornerNodes);
                }
                else analyzer.DDLogger = new DomainDecompositionLogger(subdomainPlotDirectory);
            }

            analyzer.Initialize();
            analyzer.Analyze();

            // Write crack path
            Console.WriteLine("Left Crack path:");
            foreach (var point in benchmark.LeftCrack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine("\n");
            Console.WriteLine("Right Crack path:");
            foreach (var point in benchmark.RightCrack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();

            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }
    }
}
