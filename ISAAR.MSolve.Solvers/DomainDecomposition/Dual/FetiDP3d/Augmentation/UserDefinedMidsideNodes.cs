using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public class UserDefinedMidsideNodes : IMidsideNodesSelection
    {
        private readonly Dictionary<ISubdomain, HashSet<INode>> midsideNodesOfSubdomains;
        private readonly List<INode> midsideNodesGlobal;

        public UserDefinedMidsideNodes(Dictionary<ISubdomain, HashSet<INode>> midsideNodesOfSubdomains, IDofType[] dofsPerNode)
        {
            this.midsideNodesOfSubdomains = midsideNodesOfSubdomains;
            this.DofsPerNode = dofsPerNode;
            var globalNodes = new SortedSet<INode>(); // I sort them only to match the order of Qr columns in the tests. 
            foreach (IEnumerable<INode> subdomainNodes in midsideNodesOfSubdomains.Values)
            {
                globalNodes.UnionWith(subdomainNodes);
            }
            midsideNodesGlobal = globalNodes.ToList();
        }

        public List<INode> MidsideNodesGlobal => midsideNodesGlobal;

        public IDofType[] DofsPerNode { get; }

        public HashSet<INode> GetMidsideNodesOfSubdomain(ISubdomain subdomain) => midsideNodesOfSubdomains[subdomain];

        public void Update() { }
    }
}
