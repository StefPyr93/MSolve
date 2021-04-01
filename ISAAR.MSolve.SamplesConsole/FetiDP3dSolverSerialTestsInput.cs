using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.LinearSystems;
using Xunit;
using System.Linq;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

//TODO: Perhaps I should also check intermediate steps by pulling the solver's compenent using reflection and check their state
//      and operations.
namespace ISAAR.MSolve.SamplesConsole
{
    public static class FetiDP3dSolverSerialTestsInput
    {
        //[Fact]
        public static void TestSolutionGlobalDisplacements()
        {
            (IModel model, FetiDP3dSolverSerial solver) = CreateModelAndSolver();
            RunAnalysis(model, (ISolverMpi)solver);
            Vector globalUFeti = solver.GatherGlobalDisplacements();

            Model modelDirect = ModelCreatorInput.CreateModel();
            Utilities.RemoveSubdomains(modelDirect);
            IVectorView globalUDirect = SolveDirect(modelDirect);

            Vector globalU1_2 = ReorderDirectSolverSolutionIn_globalU1_format(globalUFeti, (Vector)globalUDirect, (Model)model, modelDirect);
            Vector check = globalUFeti - globalU1_2;
            // Check solution
            double normalizedError = check.Norm2() / globalUDirect.Norm2();
            

        }

        private static Vector ReorderDirectSolverSolutionIn_globalU1_format(Vector globalU1, Vector globalU2, Model model1, Model model2)
        {
            //var freeNodesIds = model1.NodesDictionary.Values.Where(x => x.Constraints.Count == 0).Select(x=>x.ID).ToList();
            var freeNodesIds = model1.GlobalDofOrdering.GlobalFreeDofs.GetRows().Select(x => x.ID);
            Vector globalU1_2 = Vector.CreateZero(globalU1.Length);
            foreach (int nodeID in freeNodesIds)
            {
                foreach (var dof in model2.GlobalDofOrdering.GlobalFreeDofs.GetDataOfRow(model2.GetNode(nodeID)).Keys)
                {
                    int model1Dof_order = model1.GlobalDofOrdering.GlobalFreeDofs[model1.GetNode(nodeID), dof];
                    int model2Dof_order = model2.GlobalDofOrdering.GlobalFreeDofs[model2.GetNode(nodeID), dof];

                    globalU1_2[model1Dof_order] = globalU2[model2Dof_order];
                }
            }

            return globalU1_2;
        }

        //[Fact]
        //public static void TestSolutionSubdomainDisplacements()
        //{
        //    (IModel model, FetiDP3dSolverSerial solver) = CreateModelAndSolver();
        //    RunAnalysis(model, solver);

        //    // Check solution
        //    foreach (ISubdomain subdomain in model.EnumerateSubdomains())
        //    {
        //        double tol = 1E-3;
        //        IVectorView ufComputed = solver.GetLinearSystem(subdomain).Solution;
        //        Vector ufExpected = Example4x4x4Quads.ExpectedSolutions.uFree(subdomain.ID);

        //        double error = ufExpected.Subtract(ufComputed).Norm2() / ufExpected.Norm2();
        //        //Assert.Equal(0.0, 5);

        //        Assert.True(ufExpected.Equals(ufComputed, tol));
        //    }
        //}

        internal static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis
            solver.OrderDofs(false);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
            solver.BuildGlobalMatrix(new ElementStructuralStiffnessProvider());
            model.ApplyLoads();
            LoadingUtilities.ApplyNodalLoads(model, solver);
            solver.Solve();
            
        }

        private static (IModel, FetiDP3dSolverSerial) CreateModelAndSolver()
        {
            // Prepare solver
            IModel model = ModelCreatorInput.CreateModel();
            model.ConnectDataStructures();
            ICornerNodeSelection cornerNodes =ModelCreatorInput.DefineCornerNodeSelectionSerial((Model)model);
            IMidsideNodesSelection midsideNodesSelection = ModelCreatorInput.DefineMidsideNodeSelectionSerial(model);

            #region debug
            string path = @"C:\Users\Serafeim\Desktop\FETI-DP\Plots";
            var logger = new MSolve.Logging.DomainDecomposition.DomainDecompositionLoggerFetiDP(path, cornerNodes, midsideNodesSelection, true);
            //logger.PlotSubdomains(model);
            #endregion

            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = 1E-4,
                MaxIterationsProvider = new FixedMaxIterationsProvider(1000)
            };
            var matrixManagerFactory = new FetiDP3dMatrixManagerFactoryDense();
            var solverBuilder = new FetiDP3dSolverSerial.Builder(matrixManagerFactory);
            solverBuilder.PcgSettings = pcgSettings;
            solverBuilder.StiffnessDistribution = StiffnessDistributionType.HeterogeneousCondensed;
            solverBuilder.Preconditioning = new DirichletPreconditioning();
            solverBuilder.CrosspointStrategy = new FullyRedundantConstraints();

            FetiDP3dSolverSerial solver = solverBuilder.Build(model, cornerNodes,midsideNodesSelection);

            return (model, solver);
        }

        private static IVectorView SolveDirect(IModel model)
        {
            // Solver
            MSolve.Solvers.Direct.SkylineSolver solver = (new MSolve.Solvers.Direct.SkylineSolver.Builder()).BuildSolver(model);

            // Structural problem provider
            var provider = new Problems.ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new Analyzers.LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new Analyzers.StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[0].Solution;
        }
    }
}
