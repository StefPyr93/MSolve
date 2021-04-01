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
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

//TODO: Mock all other FETI classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dPreconditionerSerialTests
    {
        [Fact]
        public static void TestDiagonalDirichletPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new DiagonalDirichletPreconditioning());
            double tol = 1E-4;
            Assert.True(Example4x4x4Quads.ExpectedGlobalMatrices.PreconditionerDiagonalDirichlet.Equals(preconditioner, tol));
        }

        [Fact]
        public static void TestDirichletPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new DirichletPreconditioning());
            double tol = 1E-4;
            Assert.True(Example4x4x4Quads.ExpectedGlobalMatrices.PreconditionerDirichlet.Equals(preconditioner, tol));
        }

        [Fact]
        public static void TestLumpedPreconditioner()
        {
            Matrix preconditioner = CalcPreconditioner(new LumpedPreconditioning());
            double tol = 1E-5;
            Assert.True(Example4x4x4Quads.ExpectedGlobalMatrices.PreconditionerLumped.Equals(preconditioner, tol));
        }

        private static Matrix CalcPreconditioner(IFetiPreconditioningOperations preconditioning)
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

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
