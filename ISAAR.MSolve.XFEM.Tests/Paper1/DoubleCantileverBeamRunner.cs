using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;

namespace ISAAR.MSolve.XFEM.Tests.Paper1
{
    public class DoubleCantileverBeamRunner
    {
        private const double tipEnrichementRadius = 0.0;
        private const string crackPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\DCB\Plots\LSM";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\DCB\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\Paper1\DCB\solver_log.txt";

        public static void RunTest()
        {
            int numElementsY = 15;
            int numSubdomainsY = 3;

            // Feti-DP serial
            int numSubdomainsX = 3 * numSubdomainsY;
            bool reanalysis = false;
            bool plotSubdomains = true;
            bool plotLSM = false;
            SolveFetiDPSerial(CreateBenchmark(numElementsY, numSubdomainsX, numSubdomainsY, tipEnrichementRadius, plotLSM),
                    plotSubdomains, reanalysis);
        }

        public static void RunVariousMeshes()
        {
            // (numElementsY, dofs): (40, 10025), (75, 34664), (100, 61009), (125, 95258), (130, 102711), (200, 241900), 
            // (250, 377498), (300, 542590), (400, 963700) 
            //  (410, 1012128), (500, 1.5E6), (585, 2059800), (600, 2.165E6), (1000, 6.5E6)
            int[] numElementsY = { 75/*, 125, 250, 400, 600*/ };

            // Direct
            bool suiteSparse = true;
            bool plotLSM = true;
            //foreach (int nely in numElementsY)
            //{
            //    SolveDirect(CreateBenchmark(nely, 1, 1, tipEnrichementRadius, plotLSM), suiteSparse);
            //}

            // Feti-DP serial
            int numSubdomainsY = 25;
            int numSubdomainsX = 3 * numSubdomainsY;
            bool reanalysis = false;
            bool plotSubdomains = true;
            foreach (int nely in numElementsY)
            {
                SolveFetiDPSerial(CreateBenchmark(nely, numSubdomainsX, numSubdomainsY, tipEnrichementRadius, plotLSM),
                    plotSubdomains, reanalysis);
            }

            // Feti-DP serial reanalysis
            //reanalysis = true;
            //foreach (int nely in numElementsY)
            //{
            //    SolverFetiDPSerial(CreateBenchmark(nely, numSubdomainsX, numSubdomainsY, tipEnrichementRadius, plotLSM), 
            //        false, reanalysis);
            //}

            Console.Write("\nEnd");

        }

        public static void RunVariousSubdomains()
        {
            // SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS 
            //Make sure it performs reanalysis. I think that the dof separation is performed, even in unmodified subdomains
            // SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS 
            
            int numElementsY = 600;
            int[] numSubdomainsY = { 3, 5, 9, 15, 25 };
            bool reanalysis = true;
            bool plotLSM = false;
            bool plotSubdomains = false;
            foreach (int nsuby in numSubdomainsY)
            {
                SolveFetiDPSerial(CreateBenchmark(numElementsY, 3 * nsuby, nsuby, tipEnrichementRadius, plotLSM),
                  plotSubdomains, reanalysis);
            }
        }

        public static void RunVariousPreconditioners()
        {

        }

        public static void RunVariousProcesses()
        {

        }

        private static DcbBenchmarkBelytschko CreateBenchmark(int numElementsY, int numSubdomainsX, int numSubdomainsY, 
            double tipEnrichmentRadius, bool plotLSM)
        {
            var builder = new DcbBenchmarkBelytschko.Builder(numElementsY, numSubdomainsX, numSubdomainsY);
            builder.LsmPlotDirectory = plotLSM ? crackPlotDirectory : null;
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

        //private static ISolverMpi DefineSolver(DcbBenchmarkBelytschko benchmark, SolverType solverType)
        //{
        //    if (solverType == SolverType.FetiDP)
        //    {
        //        //benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
        //        benchmark.Model.ConnectDataStructures();

        //        //Dictionary<ISubdomain, HashSet<INode>> initialCorners = FindCornerNodesFromRectangleCorners(benchmark.Model);
        //        Func<ISubdomain, HashSet<INode>> getInitialCorners = sub => new HashSet<INode>(
        //            CornerNodeUtilities.FindCornersOfRectangle2D(sub).Where(node => node.Constraints.Count == 0));
        //        var cornerNodeSelection = new CrackedFetiDPCornerNodesSerial(benchmark.Model, benchmark.Crack, getInitialCorners);
        //        //var reordering = new OrderingAmdCSparseNet();  // This is slower than natural ordering
        //        IReorderingAlgorithm reordering = null;
        //        var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
        //        var builder = new FetiDPSolverSerial.Builder(fetiMatrices);
        //        //builder.Preconditioning = new LumpedPreconditioning();
        //        //builder.Preconditioning = new DiagonalDirichletPreconditioning();
        //        builder.Preconditioning = new DirichletPreconditioning();
        //        builder.ProblemIsHomogeneous = true;
        //        builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
        //        return builder.Build(benchmark.Model, cornerNodeSelection);
        //    }
        //    else throw new ArgumentException("Invalid solver choice.");
        //}

        private static void SolveDirect(DcbBenchmarkBelytschko benchmark, bool suiteSparse)
        {
            ISolver solver;

            if (suiteSparse)
            {
                // Suitesparse solver
                LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
                var builder = new SuiteSparseSolver.Builder();
                solver = builder.BuildSolver(benchmark.Model);
            }
            else
            {
                // Skyline solver
                var builder = new SkylineSolver.Builder();
                solver = builder.BuildSolver(benchmark.Model);

            }

            // Run analysis
            var analyzer = new QuasiStaticCrackPropagationAnalyzerOLD(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations);
            analyzer.Initialize();
            analyzer.Analyze();

            // Output
            WriteCrackPath(benchmark);
            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }

        private static void SolveFetiDPSerial(DcbBenchmarkBelytschko benchmark, bool plotSubdomains, bool reanalysis)
        {
            // Choose corner nodes
            //Dictionary<ISubdomain, HashSet<INode>> initialCorners = FindCornerNodesFromRectangleCorners(benchmark.Model);
            //Func<ISubdomain, HashSet<INode>> getInitialCorners = sub => new HashSet<INode>(
            //    CornerNodeUtilities.FindCornersOfRectangle2D(sub).Where(node => node.Constraints.Count == 0));

            // All nodes on the whole domain's boundary + crosspoints
            Func<ISubdomain, HashSet<INode>> getInitialCorners = sub => new HashSet<INode>(
                sub.EnumerateNodes().Where(n => (n.Multiplicity > 2) || benchmark.NodeIsOnBoundary(n)).
                Where(node => node.Constraints.Count == 0));
            var cornerNodeSelection = new CrackedFetiDPCornerNodesSerial(benchmark.Model, benchmark.Crack, getInitialCorners);
            benchmark.Model.ConnectDataStructures();

            // Define solver
            //var reordering = new OrderingAmdCSparseNet();  // This is slower than natural ordering
            IReorderingAlgorithm reordering = null;
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            var builder = new FetiDPSolverSerial.Builder(fetiMatrices);
            //builder.Preconditioning = new LumpedPreconditioning();
            //builder.Preconditioning = new DiagonalDirichletPreconditioning();
            builder.Preconditioning = new DirichletPreconditioning();
            builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
            FetiDPSolverSerial solver = builder.Build(benchmark.Model, cornerNodeSelection);

            // Logging
            IDomainDecompositionLogger ddLogger = null;
            bool shuffleSubdomainColors = true;
            if (plotSubdomains) ddLogger = new DomainDecompositionLoggerFetiDP(subdomainPlotDirectory, cornerNodeSelection,
                null, shuffleSubdomainColors);

            // Run analysis
            TipAdaptivePartitioner partitioner = null;
            partitioner = new TipAdaptivePartitioner(benchmark.Crack);
            var analyzer = new QuasiStaticCrackPropagationAnalyzerSerial(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, reanalysis, partitioner);
            analyzer.DDLogger = ddLogger;
            analyzer.Initialize();
            analyzer.Analyze();

            // Output
            WriteCrackPath(benchmark);
            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }

    //private static void RunCrackPropagationAnalysis(DcbBenchmarkBelytschko benchmark, ISolverMpi solver)
    //    {
    //        TipAdaptivePartitioner partitioner = null;
    //        partitioner = new TipAdaptivePartitioner(benchmark.Crack);
    //        var analyzer = new QuasiStaticCrackPropagationAnalyzerSerial(benchmark.Model, solver, benchmark.Crack, 
    //            benchmark.FractureToughness, benchmark.MaxIterations, partitioner);

    //        // Subdomain plots
    //        if (subdomainPlotDirectory != null)
    //        {
    //            if (solver is FetiDPSolverSerial fetiDP)
    //            {
    //                analyzer.DDLogger = new DomainDecompositionLoggerFetiDP(fetiDP.CornerNodes, subdomainPlotDirectory, true);
    //            }
    //            else analyzer.DDLogger = new DomainDecompositionLogger(subdomainPlotDirectory);
    //        }

    //        analyzer.Initialize();
    //        analyzer.Analyze();

    //        // Write crack path
    //        Console.WriteLine("Crack path:");
    //        foreach (var point in benchmark.Crack.CrackPath)
    //        {
    //            Console.WriteLine($"{point.X} {point.Y}");
    //        }
    //        Console.WriteLine();

    //        solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
    //        solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
    //    }

        private static void WriteCrackPath(DcbBenchmarkBelytschko benchmark)
        {
            // Write crack path
            Console.WriteLine("\nCrack path:");
            foreach (var point in benchmark.Crack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();
        }
    }
}
