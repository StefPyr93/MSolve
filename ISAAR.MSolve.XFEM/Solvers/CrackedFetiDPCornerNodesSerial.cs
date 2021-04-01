using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.XFEM.CrackGeometry;

//TODO: It is possible that some previous corner nodes become internal due to the TipAdaptivePartitioner. How to handle this?
namespace ISAAR.MSolve.XFEM.Solvers
{
    public class CrackedFetiDPCornerNodesSerial : CrackedFetiDPCornerNodesBase
    {
        private readonly IModel model;

        public CrackedFetiDPCornerNodesSerial(IModel model, ICrackDescription crack, 
            Func<ISubdomain, HashSet<INode>> getInitialCornerNodes) : base (crack, getInitialCornerNodes)
        {
            this.model = model;
        }

        public bool AreGlobalCornerNodesModified => areGlobalCornerNodesModified;

        public override HashSet<INode> GlobalCornerNodes => cornerNodesGlobal;

        public bool AreCornerNodesOfSubdomainModified(ISubdomain subdomain) => areCornerNodesModified[subdomain];

        public override HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain) => cornerNodesOfSubdomains[subdomain];

        public override void Update()
        {
            base.UpdateSubdomainsCorners(model.EnumerateSubdomains());
            base.GatherGlobalCornerNodes();
            isFirstAnalysis = false;
        }
    }
}
