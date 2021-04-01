using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public interface IHeterogeneousDistributionLoadScaling
    {
        double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load, 
            Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses);
    }
}
