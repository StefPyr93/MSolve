using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public class UsedDefinedCornerNodes : ICornerNodeSelection
    {
        private readonly Dictionary<ISubdomain, HashSet<INode>> cornerNodesOfSubdomains;
        private readonly HashSet<INode> cornerNodesGlobal;

        public UsedDefinedCornerNodes(Dictionary<ISubdomain, HashSet<INode>> cornerNodesOfSubdomains)
        {
            this.cornerNodesOfSubdomains = cornerNodesOfSubdomains;
            this.cornerNodesGlobal = new HashSet<INode>();
            foreach (IEnumerable<INode> subdomainNodes in cornerNodesOfSubdomains.Values)
            {
                cornerNodesGlobal.UnionWith(subdomainNodes);
            }
        }

        public HashSet<INode> GlobalCornerNodes => cornerNodesGlobal;

        public HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain) => cornerNodesOfSubdomains[subdomain];

        public void Update() { }
    }
}
