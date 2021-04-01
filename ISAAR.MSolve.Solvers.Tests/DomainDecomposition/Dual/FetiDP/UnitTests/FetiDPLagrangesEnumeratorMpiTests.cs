using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using Xunit;

//TODO: Mock all other classes.
//TODO: The create up to ... methods should be replaced by a class which will expose methods to calculate stuff client needs and
//      to get that stuff.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPLagrangesEnumeratorMpiTests
    {
        public static void TestBooleanMappingMatrices(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, IFetiDPDofSeparator dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) = CreateModelDofSeparatorLagrangesEnumerator(numProcesses);

            // Check
            Assert.Equal(8, lagrangesEnumerator.NumLagrangeMultipliers);
            double tolerance = 1E-13;
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                // Calculate Bpbr matrix
                ISubdomain subdomain = model.GetSubdomain(s);
                Matrix Br = lagrangesEnumerator.GetBooleanMatrix(subdomain).CopyToFullMatrix(false);
                Matrix expectedBr = Example4x4QuadsHomogeneous.GetMatrixBr(subdomain.ID);
                Assert.True(expectedBr.Equals(Br, tolerance));
            }
        }

        internal static (ProcessDistribution, IModel, FetiDPDofSeparatorMpi, LagrangeMultipliersEnumeratorMpi) 
            CreateModelDofSeparatorLagrangesEnumerator(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) =
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            // Enumerate lagranges and calculate the boolean matrices
            var crosspointStrategy = new FullyRedundantConstraints();
            var lagrangeEnumerator = new LagrangeMultipliersEnumeratorMpi(procs, model, crosspointStrategy, dofSeparator);
            lagrangeEnumerator.CalcBooleanMatrices(dofSeparator.GetRemainderDofOrdering);

            return (procs, model, dofSeparator, lagrangeEnumerator);
        }
    }
}
