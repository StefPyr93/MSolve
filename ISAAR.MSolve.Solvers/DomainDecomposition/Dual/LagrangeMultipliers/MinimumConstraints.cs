using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class MinimumConstraints : ICrosspointStrategy
    {
        public (ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus) FindSubdomainCombinations(ISubdomain[] nodeSubdomains)
        {
            int nodeMultiplicity = nodeSubdomains.Length;
            Debug.Assert(nodeMultiplicity > 1);
            int numNodeCombos = nodeSubdomains.Length - 1;
            var subdomainsPlus = new ISubdomain[numNodeCombos];
            var subdomainsMinus = new ISubdomain[numNodeCombos];
            for (int i = 0; i < numNodeCombos; ++i)
            {
                subdomainsPlus[i] = nodeSubdomains[i];
                subdomainsMinus[i] = nodeSubdomains[i + 1];
            }
            return (subdomainsPlus, subdomainsMinus);
        }
    }
}
