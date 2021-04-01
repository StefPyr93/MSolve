using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class GlobalFreeDofOrderingGeneral : IGlobalFreeDofOrdering
    {
        private readonly DofTable globalFreeDofs;
        private readonly int numGlobalFreeDofs;
        private readonly Dictionary<int, ISubdomainFreeDofOrdering> subdomainDofOrderings;

        private Dictionary<int, int[]> subdomainToGlobalDofMaps;

        public GlobalFreeDofOrderingGeneral(IModel model, int numGlobalFreeDofs, DofTable globalFreeDofs) 
        {
            this.numGlobalFreeDofs = numGlobalFreeDofs;
            this.globalFreeDofs = globalFreeDofs;

            subdomainDofOrderings = new Dictionary<int, ISubdomainFreeDofOrdering>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainDofOrderings[subdomain.ID] = subdomain.FreeDofOrdering;
            }
        }

        public DofTable GlobalFreeDofs => globalFreeDofs;

        public int NumGlobalFreeDofs => numGlobalFreeDofs;

        public void AddVectorSubdomainToGlobal(ISubdomain subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            CreateSubdomainGlobalMaps();
            int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain.ID];
            globalVector.AddIntoThisNonContiguouslyFrom(subdomainToGlobalDofs, subdomainVector);
        }

        public void AddVectorSubdomainToGlobalMeanValue(ISubdomain subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            throw new NotImplementedException();
        }

        public void ExtractVectorSubdomainFromGlobal(ISubdomain subdomain, IVectorView globalVector, IVector subdomainVector)
        {
            CreateSubdomainGlobalMaps();
            int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain.ID];
            subdomainVector.CopyNonContiguouslyFrom(globalVector, subdomainToGlobalDofs);
        }

        public ISubdomainFreeDofOrdering GetSubdomainDofOrdering(ISubdomain subdomain) => subdomainDofOrderings[subdomain.ID];

        public int[] MapSubdomainToGlobalDofs(ISubdomain subdomain)
        {
            CreateSubdomainGlobalMaps();
            return subdomainToGlobalDofMaps[subdomain.ID];
        }

        private void CreateSubdomainGlobalMaps()
        {
            if (subdomainToGlobalDofMaps == null)
            {
                subdomainToGlobalDofMaps = 
                    GlobalFreeDofOrderingUtilities.CalcSubdomainGlobalMappings(globalFreeDofs, subdomainDofOrderings);
            }
        }
    }
}
