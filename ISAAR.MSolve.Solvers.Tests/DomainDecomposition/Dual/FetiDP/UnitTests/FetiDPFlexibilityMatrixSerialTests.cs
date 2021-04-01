using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPFlexibilityMatrixSerialTests
    {
        [Fact]
        public static void TestFIrcTimesVector()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrc = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numCornerDofs, 
                flexibility.MultiplyFIrc);

            // Check
            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Equals(FIrc, tol));
        }

        [Fact]
        public static void TestFIrcTransposedTimesVector()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrcTranposed = ImplicitMatrixUtilities.MultiplyWithIdentity(numCornerDofs, numLagranges, 
                flexibility.MultiplyFIrcTransposed);

            // Check
            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Transpose().Equals(FIrcTranposed, tol));
        }

        [Fact]
        public static void TestFIrrTimesVector()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numLagranges, 
                flexibility.MultiplyFIrr);

            // Check
            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrr.Equals(FIrr, tol));
        }

        [Fact]
        public static void TestFIrrAndFIrcTransposedTimesVector()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numLagranges, 
                x => flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(x).FIrrTimesVector);
            Matrix FIrcTranposed = ImplicitMatrixUtilities.MultiplyWithIdentity(numCornerDofs, numLagranges, 
                x => flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(x).FIrcTransposedTimesVector);

            // Check
            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrr.Equals(FIrr, tol));
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Transpose().Equals(FIrcTranposed, tol));
        }
    }
}
