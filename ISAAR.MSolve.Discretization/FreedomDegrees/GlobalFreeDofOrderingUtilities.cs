using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    internal static class GlobalFreeDofOrderingUtilities
    {
        internal static Dictionary<int, int[]> CalcSubdomainGlobalMappings(DofTable globalFreeDofs,
            Dictionary<int, ISubdomainFreeDofOrdering> subdomainDofOrderings)
        {
            var subdomainToGlobalDofMaps = new Dictionary<int, int[]>(subdomainDofOrderings.Count);
            foreach (var subdomainOrderingPair in subdomainDofOrderings)
            {
                var subdomainToGlobalDofs = new int[subdomainOrderingPair.Value.NumFreeDofs];
                foreach ((INode node, IDofType dofType, int subdomainDofIdx) in subdomainOrderingPair.Value.FreeDofs)
                {
                    subdomainToGlobalDofs[subdomainDofIdx] = globalFreeDofs[node, dofType];
                }
                subdomainToGlobalDofMaps.Add(subdomainOrderingPair.Key, subdomainToGlobalDofs);
            }
            return subdomainToGlobalDofMaps;
        }
    }
}
