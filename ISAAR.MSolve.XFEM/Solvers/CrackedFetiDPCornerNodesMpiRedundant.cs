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
    /// Assumes all model data are stored in all processes. Each process will deal with the corner nodes of all subdomains. 
    /// Master process will also deal with global corner nodes.
    /// </summary>
    public class CrackedFetiDPCornerNodesMpiRedundant: CrackedFetiDPCornerNodesBase
    {
        private const int cornerNodesTag = 0;

        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public CrackedFetiDPCornerNodesMpiRedundant(ProcessDistribution processDistribution, IModel model, 
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
            else
            {
                int[] subdomainIDs = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                
                //TODO: unfortunately this does not work, due to accessing subdomains of boundary nodes.
                //IEnumerable<ISubdomain> subdomainsToUpdate = subdomainIDs.Select(s => model.GetSubdomain(s));
                IEnumerable<ISubdomain> subdomainsToUpdate = model.EnumerateSubdomains();

                base.UpdateSubdomainsCorners(subdomainsToUpdate);
            }
            isFirstAnalysis = false;
            //WriteCornerNodes();
        }

        private void WriteCornerNodes()
        {
            MpiUtilities.DoInTurn(procs.Communicator, () =>
            {
                //string path = @"C:\Users\Serafeim\Desktop\MPI\Tests\corner_nodes_process" + procs.OwnRank + ".txt";
                string path = @"C:\Users\Serafeim\Desktop\MPI\Tests\corner_nodes_parallel.txt";
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    ISubdomain sub = model.GetSubdomain(s);
                    WriteCornerNodes(sub, path, true);
                }
            });
        }
    }
}
