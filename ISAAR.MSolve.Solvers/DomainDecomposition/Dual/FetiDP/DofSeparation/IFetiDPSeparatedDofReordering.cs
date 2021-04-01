using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public interface IFetiDPSeparatedDofReordering
    {
        DofPermutation ReorderGlobalCornerDofs();
        DofPermutation ReorderSubdomainInternalDofs(ISubdomain subdomain);
        DofPermutation ReorderSubdomainRemainderDofs(ISubdomain subdomain);
    }
}