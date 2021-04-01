using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using static ISAAR.MSolve.XFEM.Tests.COMPDYN2019.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.COMPDYN2019
{
    public class DoubleCantileverBeam
    {
        private const int numElementsY = 100;
        private const double tipEnrichementRadius = 0.0;
        private const string crackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\Plots\LSM";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\solver_log.txt";

        public static void Run()
        {
            //int numSubdomainsX = 1;
            //int numSubdomainsY = 1;
            //var solverType = SolverType.Skyline;

            int numSubdomainsY = 5;
            int numSubdomainsX = 3 * numSubdomainsY;
            var solverType = SolverType.FetiDP;

            DcbBenchmarkBelytschko benchmark = CreateBenchmark(numElementsY, numSubdomainsX, numSubdomainsY, tipEnrichementRadius);
            ISolver solver = DefineSolver(benchmark, solverType);
            RunCrackPropagationAnalysis(benchmark, solver);

            Console.Write("\nEnd");
        }

        private static DcbBenchmarkBelytschko CreateBenchmark(int numElementsY, int numSubdomainsX, int numSubdomainsY, double tipEnrichmentRadius)
        {
            var builder = new DcbBenchmarkBelytschko.Builder(numElementsY, numSubdomainsX, numSubdomainsY);
            builder.LsmPlotDirectory = crackPlotDirectory;
            builder.HeavisideEnrichmentTolerance = 0.001;
            builder.MaxIterations = 8;
            builder.TipEnrichmentRadius = tipEnrichmentRadius;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            DcbBenchmarkBelytschko benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static ISolver DefineSolver(DcbBenchmarkBelytschko benchmark, SolverType solverType)
        {
            if (solverType == SolverType.Skyline)
            {
                var builder = new SkylineSolver.Builder();
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.Feti1)
            {
                benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
                double tol = 1E-2; // For denser than 50x50
                //double tol = 1E-3; // For 50x150 or coarser
                var factorizationTolerances = new Dictionary<int, double>();
                foreach (int s in benchmark.Model.Subdomains.Keys) factorizationTolerances[s] = tol;
                var fetiMatrices = new SkylineFeti1SubdomainMatrixManager.Factory();
                var builder = new Feti1Solver.Builder(fetiMatrices, factorizationTolerances);
                //builder.PreconditionerFactory = new LumpedPreconditionerOLD.Factory();
                //builder.PreconditionerFactory = new DiagonalDirichletPreconditionerOLD.Factory();
                builder.PreconditionerFactory = new DirichletPreconditionerOLD.Factory();
                builder.ProblemIsHomogeneous = true;
                var interfaceProblemSolverBuilder = new Feti1ProjectedInterfaceProblemSolver.Builder();
                interfaceProblemSolverBuilder.PcgConvergenceTolerance = 1E-7;
                builder.InterfaceProblemSolver = interfaceProblemSolverBuilder.Build();
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.FetiDP)
            {
                //benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
                benchmark.Model.ConnectDataStructures();

                Dictionary<ISubdomain, HashSet<INode>> initialCorners = FindCornerNodesFromRectangleCorners(benchmark.Model);
                var cornerNodeSelection = new CrackedFetiDPCornerNodesSerial(benchmark.Model, benchmark.Crack, 
                    sub => initialCorners[sub]);
                var fetiMatrices = new FetiDPSubdomainMatrixManagerSkylineOLD.Factory();
                var builder = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices);
                //builder.PreconditionerFactory = new LumpedPreconditionerOLD.Factory();
                //builder.PreconditionerFactory = new DiagonalDirichletPreconditionerOLD.Factory();
                builder.PreconditionerFactory = new DirichletPreconditionerOLD.Factory();
                builder.ProblemIsHomogeneous = true;
                var interfaceProblemSolverBuilder = new FetiDPInterfaceProblemSolverOLD.Builder();
                interfaceProblemSolverBuilder.PcgConvergenceTolerance = 1E-7;
                builder.InterfaceProblemSolver = interfaceProblemSolverBuilder.Build();
                return builder.BuildSolver(benchmark.Model);
            }
            else throw new ArgumentException("Invalid solver choice.");
        }

        private static void PlotSubdomains(DcbBenchmarkBelytschko benchmark, ISolver solver)
        {
            string subdomainPlotPath = subdomainPlotDirectory + "\\subdomains.vtk";
            string boundaryNodesPlotPath = subdomainPlotDirectory + "\\boundary_nodes.vtk";
            string cornerNodesPlotPath = subdomainPlotDirectory + "\\corner_nodes.vtk";

            if (solver is Feti1Solver feti1)
            {
                benchmark.Model.ConnectDataStructures();
                var writer = new MeshPartitionWriter();
                writer.WriteSubdomainElements(subdomainPlotPath, benchmark.Model);
                writer.WriteBoundaryNodes(boundaryNodesPlotPath, benchmark.Model);
            }
            else if (solver is FetiDPSolverOLD fetiDP)
            {
                benchmark.Model.ConnectDataStructures();
                var writer = new MeshPartitionWriter();
                writer.WriteSubdomainElements(subdomainPlotPath, benchmark.Model);
                writer.WriteBoundaryNodes(boundaryNodesPlotPath, benchmark.Model);

                var allCornerNodes = new HashSet<INode>();
                foreach (IEnumerable<INode> cornerNodes in fetiDP.CornerNodesOfSubdomains.Values)
                {
                    allCornerNodes.UnionWith(cornerNodes);
                }
                writer.WriteSpecialNodes(cornerNodesPlotPath, "corner_nodes", allCornerNodes);
            }
            else throw new ArgumentException("Invalid solver");
        }

        private static void RunCrackPropagationAnalysis(DcbBenchmarkBelytschko benchmark, ISolver solver)
        {
            var analyzer = new QuasiStaticCrackPropagationAnalyzerOLD(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, benchmark.Partitioner);

            // Subdomain plots
            if (subdomainPlotDirectory != null)
            {
                if (solver is FetiDPSolverOLD fetiDP)
                {
                    analyzer.DDLogger = new DomainDecompositionLoggerFetiDP(subdomainPlotDirectory, fetiDP.CornerNodeSelection, 
                        null, true);
                }
                else analyzer.DDLogger = new DomainDecompositionLogger(subdomainPlotDirectory);
            }

            analyzer.Initialize();
            analyzer.Analyze();

            // Write crack path
            Console.WriteLine("Crack path:");
            foreach (var point in benchmark.Crack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();

            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }
    }
}
