using System;
using System.Collections.Generic;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPInterfaceProblemSerialTests
    {
        [Fact]
        public static void TestInterfaceProblemMatrix()
        {
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(Example4x4QuadsHomogeneous.CreateModel());
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();

            var interfaceMatrix = new FetiDPInterfaceProblemMatrixSerial(matrixManager, flexibility);

            // Create explicit matrix that can be checked
            Matrix A = ImplicitMatrixUtilities.MultiplyWithIdentity(
                interfaceMatrix.NumRows, interfaceMatrix.NumColumns, interfaceMatrix.Multiply);

            // Check
            double tol = 1E-9;
            Assert.True(Example4x4QuadsHomogeneous.InterfaceProblemMatrix.Equals(A, tol));
        }

        [Fact]
        public static void TestInterfaceProblemRhs()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            IFetiDPMatrixManager matrixManager = new MockMatrixManager(Example4x4QuadsHomogeneous.CreateModel());
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            Vector globalDr = Example4x4QuadsHomogeneous.VectorDr;

            var interfaceSolver = new FetiDPInterfaceProblemSolverSerial(model, new PcgSettings());
            MethodInfo method = interfaceSolver.GetType().GetMethod("CalcInterfaceProblemRhs",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector pcgRhs = (Vector)method.Invoke(interfaceSolver, new object[] { matrixManager, flexibility, globalDr });

            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.InterfaceProblemRhs.Equals(pcgRhs, tol));
        }

        [Fact]
        public static void TestInterfaceProblemSolution()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            var stiffnessDistribution = new MockHomogeneousStiffnessDistribution();
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            var precondFactory = new FetiPreconditionerSerial.Factory();
            IFetiPreconditioner preconditioner = precondFactory.CreatePreconditioner(new DirichletPreconditioning(), model,
                dofSeparator, lagrangesEnumerator, matrixManager, stiffnessDistribution);

            var pcgSettings = new PcgSettings { ConvergenceTolerance = 1E-15 };
            var interfaceSolver = new FetiDPInterfaceProblemSolverSerial(model, pcgSettings);
            Vector lagranges = interfaceSolver.SolveInterfaceProblem(matrixManager,
                lagrangesEnumerator, flexibility, preconditioner, Example4x4QuadsHomogeneous.GlobalForcesNorm, 
                new SolverLoggerSerial("Test method"));

            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers.Equals(lagranges, tol));
        }

        [Fact]
        public static void TestVectorDr()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            var interfaceSolver = new FetiDPInterfaceProblemSolverSerial(model, new PcgSettings());
            MethodInfo method = interfaceSolver.GetType().GetMethod("CalcGlobalDr",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector globalDr = (Vector)method.Invoke(interfaceSolver, new object[] { matrixManager, lagrangesEnumerator });

            double tol = 1E-13;
            Assert.True(Example4x4QuadsHomogeneous.VectorDr.Equals(globalDr, tol));
        }
    }
}
