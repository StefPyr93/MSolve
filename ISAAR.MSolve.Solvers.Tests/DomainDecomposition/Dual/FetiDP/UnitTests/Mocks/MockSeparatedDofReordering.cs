using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

// Use Moq package instead of creating classes.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks
{
    public class MockSeparatedDofReordering : IFetiDPSeparatedDofReordering
    {
        public DofPermutation ReorderGlobalCornerDofs()
            => DofPermutation.CreateNoPermutation();

        public DofPermutation ReorderSubdomainInternalDofs(ISubdomain subdomain)
        => DofPermutation.CreateNoPermutation();

        public DofPermutation ReorderSubdomainRemainderDofs(ISubdomain subdomain)
            => DofPermutation.CreateNoPermutation();
    }
}
