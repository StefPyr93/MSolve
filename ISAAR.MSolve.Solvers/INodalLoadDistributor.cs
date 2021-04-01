using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers
{
    public interface INodalLoadDistributor
    {
        double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load);
    }
}
