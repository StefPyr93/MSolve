using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: The solver should decide which subdomains will be reused. This class only provides functionality.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the unconstrained freedom degrees of each subdomain and the shole model. Also applies any reordering and other 
    /// optimizations.
    /// </summary>
    public class DofOrderer : DofOrdererBase
    {
        private readonly bool doOptimizationsIfSingleSubdomain = true; // No idea why someone would want this to be false.
        
        public DofOrderer(IFreeDofOrderingStrategy freeOrderingStrategy, IDofReorderingStrategy reorderingStrategy,
            bool cacheElementToSubdomainDofMaps = true, bool doOptimizationsIfSingleSubdomain = true):
            base(freeOrderingStrategy, reorderingStrategy, cacheElementToSubdomainDofMaps)
        {
            this.doOptimizationsIfSingleSubdomain = doOptimizationsIfSingleSubdomain;
        }

        public override void OrderFreeDofs(IModel model)
        {
            if (doOptimizationsIfSingleSubdomain && (model.NumSubdomains == 1))
            {
                // Order subdomain dofs
                ISubdomain subdomain = model.EnumerateSubdomains().First();
                ISubdomainFreeDofOrdering subdomainOrdering = OrderFreeDofs(subdomain);
                subdomain.FreeDofOrdering = subdomainOrdering;

                // Order global dofs and assosiate them with subdomain dofs
                var globalOrdering = new GlobalFreeDofOrderingSingleSubdomain(subdomain, subdomainOrdering);
                model.GlobalDofOrdering = globalOrdering;
            }
            else
            {
                // Order subdomain dofs
                var subdomainOrderings = new Dictionary<ISubdomain, ISubdomainFreeDofOrdering>(model.NumSubdomains);
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    ISubdomainFreeDofOrdering subdomainOrdering = OrderFreeDofs(subdomain);
                    subdomainOrderings.Add(subdomain, subdomainOrdering);
                    subdomain.FreeDofOrdering = subdomainOrdering;
                }

                // Order global dofs
                (int numGlobalFreeDofs, DofTable globalFreeDofs) = freeOrderingStrategy.OrderGlobalDofs(model);
                var globalOrdering = new GlobalFreeDofOrderingGeneral(model, numGlobalFreeDofs, globalFreeDofs);
                model.GlobalDofOrdering = globalOrdering;
            }
        }
    }
}
