using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public interface IMidsideNodesSelection
    {
        IDofType[] DofsPerNode { get; }
        List<INode> MidsideNodesGlobal { get; }

        HashSet<INode> GetMidsideNodesOfSubdomain(ISubdomain subdomain);
    }
}
