using System;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.Transfer;
using MPI;

namespace ISAAR.MSolve.SamplesConsole.Parallel
{
    public class DistributeModelExample
    {
        private const int master = 0;

        public static void RunParallelHardcoded(string[] args)
        {
            using (new MPI.Environment(ref args))
            {
                //Console.WriteLine($"(process {Communicator.rank}) Hello World!"); // Run this to check if MPI works correctly.

                int master = 0;
                //int[] processesToSubdomains = { 0, 1, 2, 3 };
                int[][] processesToSubdomains = { new int[] { 0 }, new int[] { 1 }, new int[] { 2 }, new int[] { 3 } };
                var procs = new ProcessDistribution(Communicator.world, master, processesToSubdomains);

                // Create the model in master process
                var model = new ModelMpiCentralized(procs, Quad4PlateTest.CreateModel);
                model.ConnectDataStructures();

                // Scatter subdomain data to each process
                model.ScatterSubdomains();
                int ownSubdomainID = procs.GetSubdomainIDsOfProcess(procs.OwnRank)[0];
                ISubdomain subdomain = model.GetSubdomain(ownSubdomainID);
                //Console.WriteLine($"(process {procs.OwnRank}) Subdomain {model.GetSubdomain(procs.OwnSubdomainID).ID}");

                // Order dofs
                var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
                subdomain.FreeDofOrdering = dofOrderer.OrderFreeDofs(subdomain);

                // Create linear systems
                ILinearSystem ls = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
                ls.Reset();
                subdomain.Forces = Vector.CreateZero(ls.Size);

                // Create the stiffness matrices
                var provider = new ElementStructuralStiffnessProvider();
                var assembler = new SkylineAssembler();
                SkylineMatrix stiffness = assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, 
                    subdomain.EnumerateElements(), provider);
                ls.Matrix = stiffness;

                // Print the trace of each stiffness matrix
                double trace = ls.Matrix.Trace();
                Console.WriteLine($"(process {procs.OwnRank}) Subdomain {subdomain.ID}: trace(stiffnessMatrix) = {trace}");
            }
        }

        //public static void RunParallelWithSolver(string[] args)
        //{
        //    using (new MPI.Environment(ref args))
        //    {
        //        //Console.WriteLine($"(process {Communicator.Rank}) Hello World!"); // Run this to check if MPI works correctly.

        //        int master = 0;
        //        var procs = new ProcessDistribution(Communicator.world, master, new int[] { 0, 1, 2, 3 });

        //        // Create the model in master process
        //        var model = new ModelMpi(procs, Quad4PlateTest.CreateModel);
        //        model.ConnectDataStructures();

        //        // Setup solvers, analyzers
        //        ICornerNodeSelection cornerNodes = null;
        //        var matrixManagers = new FetiDPMatrixManagerFactorySkyline(null);
        //        var solver = new FetiDPSolverMpi(procs, model, cornerNodes, matrixManagers, true);
        //        //var problem = new ProblemStructural(model, null);
        //        var provider = new ElementStructuralStiffnessProvider();
        //        var childAnalyzer = new LinearAnalyzerMpi(model, solver, null);
        //        var parentAnalyzer = new StaticAnalyzerMpi(model, solver, provider, childAnalyzer, master);

        //        // Start the analysis
        //        parentAnalyzer.Initialize();

        //        // Print the trace of each stiffness matrix
        //        double trace = solver.LinearSystem.Matrix.Trace();
        //        Console.WriteLine($"(process {procs.OwnRank}) Subdomain {solver.Subdomain.ID}: trace(stiffnessMatrix) = {trace}");
        //    }
        //}
    }
}
