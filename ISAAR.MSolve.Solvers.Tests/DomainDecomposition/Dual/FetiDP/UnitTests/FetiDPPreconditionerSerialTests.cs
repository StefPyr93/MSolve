using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

//TODO: Mock all other FETI classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPPreconditionerSerialTests
    {
        [Fact]
        public static void TestDiagonalDirichletPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new DiagonalDirichletPreconditioning());
            double tol = 1E-13;
            Assert.True(Example4x4QuadsHomogeneous.MatrixPreconditionerDiagonalDirichlet.Equals(preconditioner, tol));
        }

        [Fact]
        public static void TestDirichletPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new DirichletPreconditioning());
            double tol = 1E-13;
            Assert.True(Example4x4QuadsHomogeneous.MatrixPreconditionerDirichlet.Equals(preconditioner, tol));
        }

        [Fact]
        public static void TestLumpedPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new LumpedPreconditioning());
            double tol = 1E-13;
            Assert.True(Example4x4QuadsHomogeneous.MatrixPreconditionerLumped.Equals(preconditioner, tol));
        }

        private static Matrix CalcPreconditioner(IFetiPreconditioningOperations preconditioning)
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IStiffnessDistribution stiffnessDistribution = new MockHomogeneousStiffnessDistribution();
            var preconditionerFactory = new FetiPreconditionerSerial.Factory();
            IFetiPreconditioner preconditioner = preconditionerFactory.CreatePreconditioner(preconditioning,
                model, dofSeparator, lagrangesEnumerator, matrixManager, stiffnessDistribution);

            // Create explicit matrices that can be checked
            int order = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix M = ImplicitMatrixUtilities.MultiplyWithIdentity(order, order, preconditioner.SolveLinearSystem);

            return M;
        }
    }
}
