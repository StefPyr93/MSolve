using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.Paper1
{
    public static class FourPointBendingBeamRunner
    {
        private const int numElementsX = 80, numElementsY = 20;
        private const double tipEnrichementRadius = 0.0;
        private const string crackPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\FourPoint\Plots\LSM";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\Paper1\FourPoint\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\Paper1\FourPoint\solver_log.txt";

        public static void Run()
        {
            SolveDirect(CreateBenchmark(numElementsX, numElementsY, 1, 1, tipEnrichementRadius));

            //int numSubdomainsX = 12, numSubdomainsY = 3;
            //FourPointBendingBeamBenchmark benchmark = CreateBenchmark(numElementsX, numElementsY, numSubdomainsX, numSubdomainsY, 
            //    tipEnrichementRadius);
            //SolveFetiDPSerial(benchmark);

            Console.Write("\nEnd");
        }

        private static FourPointBendingBeamBenchmark CreateBenchmark(int numElementsX, int numElementsY, int numSubdomainsX, 
            int numSubdomainsY, double tipEnrichmentRadius)
        {
            var builder = new FourPointBendingBeamBenchmark.Builder(numElementsX, numElementsY, numSubdomainsX, numSubdomainsY);
            builder.LsmPlotDirectory = crackPlotDirectory;
            builder.SubdomainPlotDirectory = subdomainPlotDirectory;
            builder.HeavisideEnrichmentTolerance = 0.001;
            builder.MaxIterations = 8;
            builder.TipEnrichmentRadius = tipEnrichmentRadius;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            FourPointBendingBeamBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static void SolveDirect(FourPointBendingBeamBenchmark benchmark)
        {
            var builder = new SkylineSolver.Builder();
            SkylineSolver solver = builder.BuildSolver(benchmark.Model);

            var analyzer = new QuasiStaticCrackPropagationAnalyzerOLD(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations);
            analyzer.Initialize();
            analyzer.Analyze();

            WriteCrackPath(benchmark);
            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }

        private static void SolveFetiDPSerial(FourPointBendingBeamBenchmark benchmark)
        {
            TipAdaptivePartitioner partitioner = null;
            //partitioner = new TipAdaptivePartitioner(benchmark.Crack);

            Func<ISubdomain, HashSet<INode>> getInitialCorners = sub => new HashSet<INode>(
                CornerNodeUtilities.FindCornersOfRectangle2D(sub).Where(node => node.Constraints.Count == 0));
            var cornerNodeSelection = new CrackedFetiDPCornerNodesSerial(benchmark.Model, benchmark.Crack, getInitialCorners);
            //var reordering = new OrderingAmdCSparseNet();  // This is slower than natural ordering
            IReorderingAlgorithm reordering = null;
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            var builder = new FetiDPSolverSerial.Builder(fetiMatrices);
            //builder.Preconditioning = new LumpedPreconditioning();
            //builder.Preconditioning = new DiagonalDirichletPreconditioning();
            builder.Preconditioning = new DirichletPreconditioning();
            builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
            FetiDPSolverSerial solver = builder.Build(benchmark.Model, cornerNodeSelection);
            IDomainDecompositionLogger ddLogger = null;
            if (subdomainPlotDirectory != null) ddLogger = new DomainDecompositionLoggerFetiDP(subdomainPlotDirectory, 
                cornerNodeSelection, null, true);

            var analyzer = new QuasiStaticCrackPropagationAnalyzerSerial(benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, true, partitioner);
            analyzer.DDLogger = ddLogger;
            analyzer.Initialize();
            analyzer.Analyze();

            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }

        private static void WriteCrackPath(FourPointBendingBeamBenchmark benchmark)
        {
            // Write crack path
            Console.WriteLine("Crack path:");
            foreach (var point in benchmark.Crack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();
        }
    }
}
