using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Reordering;
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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.XFEM.Analyzers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using MPI;

namespace ISAAR.MSolve.XFEM.Tests.Paper1
{
    public class DoubleCantileverBeamRunnerMpi
    {
        private const int numElementsY = 15;
        private const double tipEnrichmentRadius = 0.0;
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\Paper1\DCB\solver_log.txt";

        public static void RunTest(string[] args)
        {
            // SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS 
            // Find out why it gives slightly different results depending on the number of processes.
            // SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS SOS 

            int numSubdomainsY = 3;
            bool reanalysis = true;

            int numProcesses = int.Parse(args[0]);
            using (new MPI.Environment(ref args))
            {
                int numSubdomainsX = 3 * numSubdomainsY;

                var procs = ProcessDistribution.CreateDistribution(numProcesses, numSubdomainsX * numSubdomainsY);
                DcbBenchmarkBelytschkoMpi benchmark = CreateBenchmark(procs, numElementsY, numSubdomainsX, numSubdomainsY,
                    tipEnrichmentRadius);
                ISolverMpi solver = DefineSolver(procs, benchmark);
                RunCrackPropagationAnalysis(procs, benchmark, solver, reanalysis);
            }
        }

        private static DcbBenchmarkBelytschkoMpi CreateBenchmark(ProcessDistribution procs, int numElementsY, int numSubdomainsX,
            int numSubdomainsY, double tipEnrichmentRadius)
        {
            var builder = new DcbBenchmarkBelytschkoMpi.Builder(procs, numElementsY, numSubdomainsX, numSubdomainsY);
            //builder.LsmPlotDirectory = crackPlotDirectory;
            //builder.SubdomainPlotDirectory = subdomainPlotDirectory;
            builder.HeavisideEnrichmentTolerance = 0.001;
            builder.MaxIterations = 8;
            builder.TipEnrichmentRadius = tipEnrichmentRadius;
            builder.JintegralRadiusOverElementSize = 2.0;

            DcbBenchmarkBelytschkoMpi benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static ISolverMpi DefineSolver(ProcessDistribution procs, DcbBenchmarkBelytschkoMpi benchmark)
        {
            benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
            Func<ISubdomain, HashSet<INode>> getInitialCorners = sub => new HashSet<INode>(
                    CornerNodeUtilities.FindCornersOfRectangle2D(sub).Where(node => node.Constraints.Count == 0));
            //var cornerNodeSelection = new CrackedFetiDPCornerNodesMpiCentralized(procs, benchmark.Model, benchmark.Crack, 
            //    getInitialCorners);
            var cornerNodeSelection = new CrackedFetiDPCornerNodesMpiRedundant(procs, benchmark.Model, benchmark.Crack,
                getInitialCorners);
            //var reordering = new OrderingAmdCSparseNet();  // This is slower than natural ordering
            IReorderingAlgorithm reordering = null;
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            var builder = new FetiDPSolverMpi.Builder(procs, fetiMatrices);
            //builder.Preconditioning = new LumpedPreconditioning();
            //builder.Preconditioning = new DiagonalDirichletPreconditioning();
            builder.Preconditioning = new DirichletPreconditioning();
            builder.ProblemIsHomogeneous = true;
            builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
            return builder.Build(benchmark.Model, cornerNodeSelection);
        }

        private static void RunCrackPropagationAnalysis(ProcessDistribution procs, DcbBenchmarkBelytschkoMpi benchmark,
            ISolverMpi solver, bool reanalysis)
        {
            TipAdaptivePartitioner partitioner = null;
            partitioner = new TipAdaptivePartitioner(benchmark.Crack);
            //var analyzer = new QuasiStaticCrackPropagationAnalyzerMpiCentralized(procs, benchmark.Model, solver, benchmark.Crack,
            //    benchmark.FractureToughness, benchmark.MaxIterations, partitioner);
            var analyzer = new QuasiStaticCrackPropagationAnalyzerMpiRedundnat(procs, benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, reanalysis, partitioner);

            analyzer.Initialize();
            analyzer.Analyze();

            if (procs.IsMasterProcess)
            {
                solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
                solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);

                // Write crack path
                Console.WriteLine("Crack path:");
                foreach (var point in benchmark.Crack.SingleCracks[0].CrackPath)
                {
                    Console.WriteLine($"{point.X} {point.Y}");
                }
                Console.WriteLine();
            }
        }
    }
}
