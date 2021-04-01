using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class SubdomainLagrangeMultiplier
    {
        internal SubdomainLagrangeMultiplier(int globalLagrangeIndex, INode node, IDofType dof, bool subdomainSign)
        {
            this.GlobalLagrangeIndex = globalLagrangeIndex;
            this.Node = node;
            this.DofType = dof;
            this.SubdomainSign = subdomainSign;
        }

        internal IDofType DofType { get; }
        internal int GlobalLagrangeIndex { get; }
        internal INode Node { get; }
        internal bool SubdomainSign { get; }
    }
}
