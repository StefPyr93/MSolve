using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Is it even necessary to have this? It does not coordinate anything. Can't I just pass the global/subdomain ones to the corresponding methods?
namespace ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation
{
    public interface IDofSeparator
    {
        /// <summary>
        /// Dofs where Lagrange multipliers will be applied. Depending on the domain decomposition method, these could be all
        /// free boundary dofs (e.g. FETI-1) or a subset of them (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; }


        /// <summary>
        /// The indices of dofs, where Lagrange multipliers will be applied, into a sequence of subdomain dofs. Depending on the 
        /// domain decomposition method, that sequence could be all free boundary dofs (e.g. FETI-1) or a subset of them 
        /// (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        int[] GetBoundaryDofIndices(ISubdomain subdomain);

        /// <summary>
        /// Dofs of each subdomain where Lagrange multipliers will be applied. Depending on the domain decomposition method, 
        /// these could be all free boundary dofs (e.g. FETI-1) or a subset of them (e.g. boundary dofs minus corner dofs in 
        /// FETI-DP).
        /// </summary>
        (INode node, IDofType dofType)[] GetBoundaryDofs(ISubdomain subdomain);


        /// <summary>
        /// The indices of dofs, which only occur for each subdomain, into a sequence of subdomain dofs. Depending on the 
        /// domain decomposition method, that sequence could be all free boundary dofs (e.g. FETI-1) or a subset of them 
        /// (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        int[] GetInternalDofIndices(ISubdomain subdomain);
    }
}
