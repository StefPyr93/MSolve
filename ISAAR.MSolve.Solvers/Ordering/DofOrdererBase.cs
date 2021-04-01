using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Ordering
{
    public abstract class DofOrdererBase : IDofOrderer
    {
        //TODO: this should also be a strategy, so that I could have caching with fallbacks, in case of insufficient memor.
        protected readonly bool cacheElementToSubdomainDofMaps = true;
        protected readonly ConstrainedDofOrderingStrategy constrainedOrderingStrategy;
        protected readonly IFreeDofOrderingStrategy freeOrderingStrategy;
        protected readonly IDofReorderingStrategy reorderingStrategy;

        public DofOrdererBase(IFreeDofOrderingStrategy freeOrderingStrategy, IDofReorderingStrategy reorderingStrategy,
            bool cacheElementToSubdomainDofMaps = true)
        {
            this.constrainedOrderingStrategy = new ConstrainedDofOrderingStrategy();
            this.freeOrderingStrategy = freeOrderingStrategy;
            this.reorderingStrategy = reorderingStrategy;
            this.cacheElementToSubdomainDofMaps = cacheElementToSubdomainDofMaps;
        }

        public ISubdomainConstrainedDofOrdering OrderConstrainedDofs(ISubdomain subdomain)
        {
            (int numConstrainedDofs, DofTable constrainedDofs) =
                constrainedOrderingStrategy.OrderSubdomainDofs(subdomain);
            if (cacheElementToSubdomainDofMaps)
            {
                return new SubdomainConstrainedDofOrderingCaching(numConstrainedDofs, constrainedDofs);
            }
            else return new SubdomainConstrainedDofOrderingGeneral(numConstrainedDofs, constrainedDofs);
        }

        public ISubdomainFreeDofOrdering OrderFreeDofs(ISubdomain subdomain)
        {
            if (!subdomain.ConnectivityModified) return subdomain.FreeDofOrdering;

            // Order subdomain dofs
            (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) = freeOrderingStrategy.OrderSubdomainDofs(subdomain);
            ISubdomainFreeDofOrdering subdomainOrdering;
            if (cacheElementToSubdomainDofMaps) subdomainOrdering = new SubdomainFreeDofOrderingCaching(
                numSubdomainFreeDofs, subdomainFreeDofs);
            else subdomainOrdering = new SubdomainFreeDofOrderingGeneral(numSubdomainFreeDofs, subdomainFreeDofs);

            // Reorder subdomain dofs
            reorderingStrategy.ReorderDofs(subdomain, subdomainOrdering);

            return subdomainOrdering;
        }

        public abstract void OrderFreeDofs(IModel model);
    }
}
