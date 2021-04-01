using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public class UsedDefinedCornerNodesMpi : ICornerNodeSelection
    {
        private readonly HashSet<INode> cornerNodesGlobal_master;
        private readonly Dictionary<ISubdomain, HashSet<INode>> cornerNodesOfSubdomains;
        private readonly ProcessDistribution procs;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="processDistribution"></param>
        /// <param name="cornerNodesOfSubdomains">
        /// For master process it must contain all subdomain data. For every other process it should only contain the data of
        /// the corresponding subdomains.
        /// </param>
        public UsedDefinedCornerNodesMpi(ProcessDistribution processDistribution, 
            Dictionary<ISubdomain, HashSet<INode>> cornerNodesOfSubdomains)
        {
            this.procs = processDistribution;
            this.cornerNodesOfSubdomains = cornerNodesOfSubdomains;
            if (processDistribution.IsMasterProcess)
            {
                this.cornerNodesGlobal_master = new HashSet<INode>();
                foreach (IEnumerable<INode> subdomainNodes in cornerNodesOfSubdomains.Values)
                {
                    cornerNodesGlobal_master.UnionWith(subdomainNodes);
                }
            }
        }

        public HashSet<INode> GlobalCornerNodes
        {
            get
            {
                procs.CheckProcessIsMaster();
                return cornerNodesGlobal_master;
            }
        }

        public HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomain.ID);
            return cornerNodesOfSubdomains[subdomain];
        }

        public void Update() { }
    }
}
