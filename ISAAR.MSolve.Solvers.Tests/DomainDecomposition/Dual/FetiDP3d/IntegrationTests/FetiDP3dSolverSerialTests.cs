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

//TODO: Perhaps I should also check intermediate steps by pulling the solver's compenent using reflection and check their state
//      and operations.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.IntegrationTests
{
    public static class FetiDP3dSolverSerialTests
    {
        [Fact]
        public static void TestSolutionGlobalDisplacements()
        {
            (IModel model, FetiDP3dSolverSerial solver) = CreateModelAndSolver();
            RunAnalysis(model, solver);
            Vector globalUFeti = solver.GatherGlobalDisplacements();

            Model modelDirect = Example4x4x4Quads.ModelCreator.CreateModel();
            Utilities.RemoveSubdomains(modelDirect);
            IVectorView globalUDirect = SolveDirect(modelDirect);

            // Check solution
            double normalizedError = globalUDirect.Subtract(globalUFeti).Norm2() / globalUDirect.Norm2();
            Assert.Equal(0.0, normalizedError, 5);

            //// Check solution
            //double tol = 1E-8;
            //Assert.True(Example4x4x4Quads.SolutionGlobalDisplacements.Equals(globalU, tol));
        }

        [Fact]
        public static void TestSolutionSubdomainDisplacements()
        {
            (IModel model, FetiDP3dSolverSerial solver) = CreateModelAndSolver();
            RunAnalysis(model, solver);

            // Check solution
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                double tol = 1E-3;
                IVectorView ufComputed = solver.GetLinearSystem(subdomain).Solution;
                Vector ufExpected = Example4x4x4Quads.ExpectedSolutions.uFree(subdomain.ID);

                double error = ufExpected.Subtract(ufComputed).Norm2() / ufExpected.Norm2();
                //Assert.Equal(0.0, 5);

                Assert.True(ufExpected.Equals(ufComputed, tol));
            }
        }

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
            IModel model = Example4x4x4Quads.ModelCreator.CreateModel();
            model.ConnectDataStructures();
            ICornerNodeSelection cornerNodes = Example4x4x4Quads.ModelCreator.DefineCornerNodeSelectionSerial(model);
            IMidsideNodesSelection midsideNodesSelection = Example4x4x4Quads.ModelCreator.DefineMidsideNodeSelectionSerial(model);

            #region debug
            string path = @"C:\Users\Serafeim\Desktop\FETI-DP\Plots";
            var logger = new MSolve.Logging.DomainDecomposition.DomainDecompositionLoggerFetiDP(path, cornerNodes, midsideNodesSelection, true);
            //logger.PlotSubdomains(model);
            #endregion

            var matrixManagerFactory = new FetiDP3dMatrixManagerFactoryDense();
            var solverBuilder = new FetiDP3dSolverSerial.Builder(matrixManagerFactory);
            FetiDP3dSolverSerial solver = solverBuilder.Build(model, cornerNodes,midsideNodesSelection);

            return (model, solver);
        }

        private static IVectorView SolveDirect(IModel model)
        {
            // Solver
            Direct.SkylineSolver solver = (new Direct.SkylineSolver.Builder()).BuildSolver(model);

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
