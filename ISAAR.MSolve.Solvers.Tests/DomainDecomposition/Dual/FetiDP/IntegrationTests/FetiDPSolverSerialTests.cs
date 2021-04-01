using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;
using Xunit;

//TODO: Perhaps I should also check intermediate steps by pulling the solver's compenent using reflection and check their state
//      and operations.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    public static class FetiDPSolverSerialTests
    {
        [Theory]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestSolutionGlobalDisplacements(MatrixFormat format)
        {
            (IModel model, FetiDPSolverSerial solver) = CreateModelAndSolver(format);
            RunAnalysis(model, solver);
            Vector globalU = solver.GatherGlobalDisplacements();

            // Check solution
            double tol = 1E-8;
            Assert.True(Example4x4QuadsHomogeneous.SolutionGlobalDisplacements.Equals(globalU, tol));
        }

        [Theory]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestSolutionSubdomainDisplacements(MatrixFormat format)
        {
            (IModel model, FetiDPSolverSerial solver) = CreateModelAndSolver(format);
            RunAnalysis(model, solver);

            // Check solution
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                double tol = 1E-6;
                IVectorView ufComputed = solver.GetLinearSystem(subdomain).Solution;
                Vector ufExpected = Example4x4QuadsHomogeneous.GetSolutionFreeDisplacements(subdomain.ID);
                Assert.True(ufExpected.Equals(ufComputed, tol));
            }
        }

        internal static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis //BMA
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

        private static (IModel, FetiDPSolverSerial) CreateModelAndSolver(MatrixFormat format)
        {
            // Prepare solver
            IModel model = Example4x4QuadsHomogeneous.CreateModel();
            model.ConnectDataStructures();
            ICornerNodeSelection cornerNodes = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            IFetiDPMatrixManagerFactory fetiMatrices = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
            FetiDPSolverSerial solver = solverBuilder.Build(model, cornerNodes);

            return (model, solver);
        }
    }
}
