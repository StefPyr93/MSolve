using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.XFEM.CrackGeometry;

//TODO: Remove duplication between this and the serial implementation.
//TODO: It is possible that some previous corner nodes become internal due to the TipAdaptivePartitioner. How to handle this?
namespace ISAAR.MSolve.XFEM.Solvers
{
    /// <summary>
    /// Assumes 1 subdomain per process. Also assumes that all model data exist in master process, while the other processes
    /// only store the entities of their corresponding subdomains.
    /// </summary>
    public class CrackedFetiDPCornerNodesMpiCentralized : CrackedFetiDPCornerNodesBase
    {
        private const int cornerNodesTag = 0;

        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public CrackedFetiDPCornerNodesMpiCentralized(ProcessDistribution processDistribution, IModel model,
            ICrackDescription crack, Func<ISubdomain, HashSet<INode>> getInitialCornerNodes) :
            base(crack, getInitialCornerNodes)
        {
            this.procs = processDistribution;
            this.model = model;
        }

        public bool AreGlobalCornerNodesModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return areGlobalCornerNodesModified;
            }
        }

        public override HashSet<INode> GlobalCornerNodes
        {
            get
            {
                procs.CheckProcessIsMaster();
                return cornerNodesGlobal;
            }
        }

        public bool AreCornerNodesOfSubdomainModified(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomain.ID);
            return areCornerNodesModified[subdomain];
        }

        public override HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomain.ID);
            return cornerNodesOfSubdomains[subdomain];
        }

        public override void Update()
        {
            if (procs.IsMasterProcess)
            {
                base.UpdateSubdomainsCorners(model.EnumerateSubdomains());
                base.GatherGlobalCornerNodes();
            }
            ScatterSubdomainCornerNodes();
            isFirstAnalysis = false;
            //WriteCornerNodes();
        }

        private void ScatterSubdomainCornerNodes()
        {
            if (procs.IsMasterProcess)
            {
                // Send to each process the corner nodes of its subdomain, if they are modified.
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) continue;
                    int s = procs.GetSubdomainIDsOfProcess(p).First();
                    ISubdomain subdomain = model.GetSubdomain(s);
                    if (subdomain.ConnectivityModified)
                    {
                        HashSet<INode> cornerNodes = cornerNodesOfSubdomains[subdomain];
                        int[] cornerIDs = cornerNodes.Select(n => n.ID).ToArray();
                        MpiUtilities.SendArray<int>(procs.Communicator, cornerIDs, p, cornerNodesTag);
                    }
                }
            }
            else
            {
                int[] subdomainIds = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                ISubdomain subdomain = model.GetSubdomain(subdomainIds.First());
                // Receive the corner nodes from master, if they are modified.
                if (subdomain.ConnectivityModified)
                {
                    int[] cornerIDs = MpiUtilities.ReceiveArray<int>(procs.Communicator, procs.MasterProcess, cornerNodesTag);
                    var cornerNodes = new HashSet<INode>();
                    foreach (int n in cornerIDs) cornerNodes.Add(subdomain.GetNode(n));
                    cornerNodesOfSubdomains[subdomain] = cornerNodes;
                }
            }
        }

        private void WriteCornerNodes()
        {
            MpiUtilities.DoInTurn(procs.Communicator, () =>
            {
                int s = procs.GetSubdomainIDsOfProcess(procs.OwnRank).First();
                Console.Write($"Process {procs.OwnRank}: Corner nodes of subdomain {s}: ");
                foreach (INode node in cornerNodesOfSubdomains[model.GetSubdomain(s)]) Console.Write(node.ID + " ");
                Console.WriteLine();
            });
        }
    }
}
