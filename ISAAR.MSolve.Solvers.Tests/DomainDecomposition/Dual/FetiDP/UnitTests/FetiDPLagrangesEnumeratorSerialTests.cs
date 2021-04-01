using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using Xunit;

//TODO: Mock all other classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPLagrangesEnumeratorSerialTests
    {
        [Fact]
        public static void TestBooleanMappingMatrices()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangeEnumerator) =
                CreateModelDofSeparatorLagrangesEnumerator();

            // Check
            Assert.Equal(8, lagrangeEnumerator.NumLagrangeMultipliers);
            double tolerance = 1E-13;
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                Matrix Br = lagrangeEnumerator.GetBooleanMatrix(subdomain).CopyToFullMatrix(false);
                Matrix expectedBr = Example4x4QuadsHomogeneous.GetMatrixBr(subdomain.ID);
                Assert.True(expectedBr.Equals(Br, tolerance));
            }
        }

        internal static (IModel, FetiDPDofSeparatorSerial, LagrangeMultipliersEnumeratorSerial)
            CreateModelDofSeparatorLagrangesEnumerator()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator) =
                FetiDPDofSeparatorSerialTests.CreateModelAndDofSeparator();

            // Enumerate lagranges and calculate the boolean matrices
            var crosspointStrategy = new FullyRedundantConstraints();
            var lagrangeEnumerator = new LagrangeMultipliersEnumeratorSerial(model, crosspointStrategy, dofSeparator);
            lagrangeEnumerator.CalcBooleanMatrices(dofSeparator.GetRemainderDofOrdering);

            return (model, dofSeparator, lagrangeEnumerator);
        }
    }
}
