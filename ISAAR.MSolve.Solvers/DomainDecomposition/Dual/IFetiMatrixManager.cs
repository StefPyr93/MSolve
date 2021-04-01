using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public interface IFetiMatrixManager
    {
        IFetiSubdomainMatrixManager GetSubdomainMatrixManager(ISubdomain subdomain);
    }
}
