using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using Xunit;

//TODO: Mock all other classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class HomogeneousStiffnessDistributionSerialTests
    {
        [Fact]
        public static void TestBooleanMappingMatrices()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            // Set up stiffness distribution
            var stiffnessDistribution = new HomogeneousStiffnessDistributionSerial(model, dofSeparator,
                new FetiDPHomogeneousDistributionLoadScaling(dofSeparator));
            stiffnessDistribution.Update();

            double tol = 1E-13;
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                // Calculate Bpbr
                SignedBooleanMatrixColMajor Bb = lagrangesEnumerator.GetBooleanMatrix(subdomain).GetColumns(
                    dofSeparator.GetBoundaryDofIndices(subdomain), false);
                IMappingMatrix Bpbr =
                    stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrix(lagrangesEnumerator, subdomain, Bb);

                // Check Bpbr
                Matrix explicitBpr = Bpbr.MultiplyRight(Matrix.CreateIdentity(Bpbr.NumColumns));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixBpbr(subdomain.ID).Equals(explicitBpr, tol));
            }
        }
    }
}
