using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using Xunit;

//TODO: Mock all other classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class HomogeneousStiffnessDistributionMpiTests
    {
        public static void TestBooleanMappingMatrices(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, IFetiDPDofSeparator dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) = 
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Caclulate scaling coefficients
            var stiffnessDistribution = new HomogeneousStiffnessDistributionMpi(procs, model, dofSeparator,
                new FetiDPHomogeneousDistributionLoadScaling(dofSeparator));
            stiffnessDistribution.Update();

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                // Calculate Bpbr matrix
                ISubdomain subdomain = model.GetSubdomain(s);
                int[] boundaryDofs = dofSeparator.GetBoundaryDofIndices(subdomain);
                SignedBooleanMatrixColMajor Bb = lagrangesEnumerator.GetBooleanMatrix(subdomain).GetColumns(boundaryDofs, false);
                IMappingMatrix Bpbr =
                    stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrix(lagrangesEnumerator, subdomain, Bb);

                // Check Bpbr matrix
                double tol = 1E-13;
                Matrix explicitBpr = Bpbr.MultiplyRight(Matrix.CreateIdentity(Bpbr.NumColumns));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixBpbr(subdomain.ID).Equals(explicitBpr, tol));
            }
        }
    }
}
