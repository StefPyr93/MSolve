using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;

//TODO: Needs a proper name. This probably cannot be incorporated in the ISubdomainDofOrdering, as the intent is different and
//      depending on the DD method the dof categories may be different (e.g. FETI-1: internal/boundary, 
//      FETI-DP: corner/boundary/remainder)
//TODO: Not sure about having the indexing data of all subdomains in the same class.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.DofSeparation
{
    public class Feti1DofSeparator : IDofSeparator
    {
        public Dictionary<int, int[]> BoundaryDofIndices { get; private set; }
        public Dictionary<int, (INode node, IDofType dofType)[]> BoundaryDofs { get; private set; }
        public Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; private set; }
        public Dictionary<int, int[]> InternalDofIndices { get; private set; }

        public Feti1DofSeparator()
        {
            InternalDofIndices = new Dictionary<int, int[]>();
            BoundaryDofIndices = new Dictionary<int, int[]>();
            BoundaryDofs = new Dictionary<int, (INode node, IDofType dofType)[]>();
        }

        public void DefineGlobalBoundaryDofs(IModel model)
        {
            GlobalBoundaryDofs = DofSeparationUtilities.DefineGlobalBoundaryDofs(model.EnumerateNodes(), 
                model.GlobalDofOrdering.GlobalFreeDofs);
        }

        public void SeparateBoundaryInternalDofs(ISubdomain subdomain)
        {
            int s = subdomain.ID;
            (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                = DofSeparationUtilities.SeparateBoundaryInternalDofs(subdomain.EnumerateNodes(), subdomain.FreeDofOrdering.FreeDofs);

            InternalDofIndices[s] = internalDofIndices;
            BoundaryDofIndices[s] = boundaryDofIndices;
            BoundaryDofs[s] = boundaryDofConnectivities;
        }

        public int[] GetBoundaryDofIndices(ISubdomain subdomain) => BoundaryDofIndices[subdomain.ID];

        public (INode node, IDofType dofType)[] GetBoundaryDofs(ISubdomain subdomain) => BoundaryDofs[subdomain.ID];

        public int[] GetInternalDofIndices(ISubdomain subdomain) => InternalDofIndices[subdomain.ID];
    }
}
