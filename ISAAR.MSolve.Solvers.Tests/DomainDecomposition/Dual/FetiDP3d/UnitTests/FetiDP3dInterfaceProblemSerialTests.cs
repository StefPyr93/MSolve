using System;
using System.Collections.Generic;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dInterfaceProblemSerialTests
    {
        [Fact]
        public static void TestInterfaceProblemMatrix()
        {
            IModel model = Example4x4x4Quads.ModelCreator.CreateModel();
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();

            var interfaceMatrix = new FetiDPInterfaceProblemMatrixSerial(matrixManager, flexibility);

            // Create explicit matrix that can be checked
            Matrix A = ImplicitMatrixUtilities.MultiplyWithIdentity(
                interfaceMatrix.NumRows, interfaceMatrix.NumColumns, interfaceMatrix.Multiply);

            // Check
            double tol = 1E-8;
            Assert.True(Example4x4x4Quads.ExpectedGlobalMatrices.InterfaceProblemMatrixSimple.Equals(A, tol));
        }

        [Fact]
        public static void TestInterfaceProblemRhs()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints = 
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            
            Vector globalDr = Example4x4x4Quads.ExpectedGlobalMatrices.VectorDr;

            var interfaceSolver = new FetiDP3dInterfaceProblemSolverSerial(model, new PcgSettings() { ConvergenceTolerance = 1E-9 }, augmentationConstraints);
            MethodInfo method = interfaceSolver.GetType().GetMethod("CalcInterfaceProblemRhs",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector pcgRhs = (Vector)method.Invoke(interfaceSolver, new object[] { matrixManager, flexibility, globalDr });

            double tol = 1E-8;
            Assert.True(Example4x4x4Quads.ExpectedGlobalMatrices.InterfaceProblemRhsSimple.Equals(pcgRhs, tol));
        }

        [Fact]
        public static void TestInterfaceProblemSolution()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            IFetiPreconditioner preconditioner = new MockPreconditioner(); //TODO: Also mock the preconditioner in 2D FETI-DP

            var pcgSettings = new PcgSettings 
            { 
                ConvergenceTolerance = Example4x4x4Quads.ExpectedGlobalMatrices.PcgResidualNormRatio, 
                MaxIterationsProvider = new FixedMaxIterationsProvider(Example4x4x4Quads.ExpectedGlobalMatrices.PcgIterations) 
            };
            var interfaceSolver = new FetiDP3dInterfaceProblemSolverSerial(model, pcgSettings, augmentationConstraints);
            Vector lagranges = interfaceSolver.SolveInterfaceProblem(matrixManager,
                lagrangesEnumerator, flexibility, preconditioner, Example4x4x4Quads.ExpectedGlobalMatrices.GlobalForcesNorm, 
                new SolverLoggerSerial("Test method"));

            double tol = 1E-8;
            Vector lagrangesExpected = Example4x4x4Quads.ExpectedSolutions.SolutionLagrangesSimple();
            Assert.True(lagrangesExpected.Equals(lagranges, tol));
        }
    }
}
