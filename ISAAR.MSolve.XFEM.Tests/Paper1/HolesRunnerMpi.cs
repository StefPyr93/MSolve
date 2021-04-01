using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Distributed;
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
    public class HolesRunnerMpi
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
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\Paper1\Holes\solver_log.txt";

        public static void RunTest(string[] args)
        {
            int numSubdomains = 10;
            bool reanalysis = true;

            int numProcesses = int.Parse(args[0]);
            using (new MPI.Environment(ref args))
            {
                var procs = ProcessDistribution.CreateDistribution(numProcesses, numSubdomains);

                // FETI-DP 10 subdomains
                HolesBenchmarkMpi benchmarkSub10 = CreateBenchmark(procs, numSubdomains);
                ISolverMpi solverFetiDP = DefineSolver(procs, benchmarkSub10);
                RunCrackPropagationAnalysis(procs, benchmarkSub10, solverFetiDP, reanalysis);

                // FETI-DP 15 subdomains
                //HolesBenchmarkMpi benchmarkSub15 = CreateMultiSubdomainBenchmark(numSubdomains);
                //ISolver solverFetiDP = DefineSolver(benchmarkSub15);
                //RunCrackPropagationAnalysis(benchmarkSub15, solverFetiDP, reanalysis);

                Console.Write("\nEnd");
            }
        }

        private static HolesBenchmarkMpi CreateBenchmark(ProcessDistribution procs, int numSubdomains)
        {
            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.

            var builder = new HolesBenchmarkMpi.Builder(procs, meshPath, growthLength);
            builder.HeavisideEnrichmentTolerance = 0.12;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // If you modify the following two parameters significantly, then you will need to redefine which nodes are expected 
            // to be enriched.
            builder.TipEnrichmentRadius = 0.5;
            builder.BC = HolesBenchmarkMpi.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;

            builder.MaxIterations = 12;

            HolesBenchmarkMpi benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel(numSubdomains);
            return benchmark;
        }

        private static ISolverMpi DefineSolver(ProcessDistribution procs, HolesBenchmarkMpi benchmark)
        {
            Func<ISubdomain, HashSet<INode>> getInitialCorners =
                (sub) => new HashSet<INode>(sub.EnumerateNodes().Where(node => node.Multiplicity > 2));
            //var cornerNodeSelection = new CrackedFetiDPCornerNodesMpiCentralized(procs, benchmark.Model, benchmark.Crack, 
            //    getInitialCorners);
            var cornerNodeSelection = new CrackedFetiDPCornerNodesMpiRedundant(procs, benchmark.Model, benchmark.Crack,
                getInitialCorners);

            IReorderingAlgorithm reordering = null;
            reordering = new OrderingAmdCSparseNet();

            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            var builder = new FetiDPSolverMpi.Builder(procs, fetiMatrices);
            //builder.Preconditioning = new LumpedPreconditioning();
            //builder.Preconditioning = new DiagonalDirichletPreconditioning();
            builder.Preconditioning = new DirichletPreconditioning();
            builder.ProblemIsHomogeneous = true;
            builder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-7 };
            return builder.Build(benchmark.Model, cornerNodeSelection);
        }

        private static void RunCrackPropagationAnalysis(ProcessDistribution procs, HolesBenchmarkMpi benchmark, 
            ISolverMpi solver, bool reanalysis)
        {
            TipAdaptivePartitioner partitioner = null;
            partitioner = new TipAdaptivePartitioner(benchmark.Crack);
            var analyzer = new QuasiStaticCrackPropagationAnalyzerMpiRedundnat(procs, benchmark.Model, solver, benchmark.Crack,
                benchmark.FractureToughness, benchmark.MaxIterations, reanalysis, partitioner);
            analyzer.Initialize();
            analyzer.Analyze();

            if (procs.IsMasterProcess)
            {
                solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
                solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);

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
            }
        }
    }
}
