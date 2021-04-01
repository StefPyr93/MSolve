using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public interface IHomogeneousDistributionLoadScaling
    {
        double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load);
    }
}
