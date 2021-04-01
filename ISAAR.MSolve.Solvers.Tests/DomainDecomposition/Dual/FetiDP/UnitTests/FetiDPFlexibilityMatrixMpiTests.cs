using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

// perhaps add tests for repeated use to simulate PCG
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPFlexibilityMatrixMpiTests
    {
        public static void TestFIrcTimesVector(int numProcesses) 
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator, 
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixMpi(procs, model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrc = ImplicitMatrixUtilities.MultiplyWithIdentityMpi(numLagranges, numCornerDofs, 
                flexibility.MultiplyFIrc);
            
            if (procs.IsMasterProcess)
            {
                // Check
                double tol = 1E-11;
                Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Equals(FIrc, tol));
            }
        }

        public static void TestFIrcTransposedTimesVector(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixMpi(procs, model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrcTransposed = ImplicitMatrixUtilities.MultiplyWithIdentityMpi(numLagranges, numCornerDofs, 
                flexibility.MultiplyFIrcTransposed);

            if (procs.IsMasterProcess)
            {
                // Check
                double tol = 1E-11;
                Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Transpose().Equals(FIrcTransposed, tol));
            }
        }

        public static void TestFIrrTimesVector(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixMpi(procs, model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentityMpi(numLagranges, numLagranges, 
                flexibility.MultiplyFIrr);

            if (procs.IsMasterProcess)
            {
                // Check
                double tol = 1E-11;
                Assert.True(Example4x4QuadsHomogeneous.MatrixFIrr.Equals(FIrr, tol));
            }
        }

        public static void TestFIrrAndFIrcTransposedTimesVector(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Setup matrix manager
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDPFlexibilityMatrixMpi(procs, model, dofSeparator, lagrangesEnumerator, matrixManager);
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentityMpi(numLagranges, numLagranges,
                x => flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(x).FIrrTimesVector);
            Matrix FIrcTransposed = ImplicitMatrixUtilities.MultiplyWithIdentityMpi(numLagranges, numCornerDofs,
                x => flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(x).FIrcTransposedTimesVector);

            if (procs.IsMasterProcess)
            {
                // Check
                double tol = 1E-11;
                Assert.True(Example4x4QuadsHomogeneous.MatrixFIrr.Equals(FIrr, tol));
                Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Transpose().Equals(FIrcTransposed, tol));
            }
        }
    }
}
