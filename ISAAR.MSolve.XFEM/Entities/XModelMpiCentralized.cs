using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Transfer;


//TODO: Transfer level sets and recalculate nodal enrichments in each process. Perhaps also identify which nodes are enriched 
//      with what. 
namespace ISAAR.MSolve.XFEM.Entities
{
    public class XModelMpiCentralized : ModelMpiCentralizedBase<XModel>, IXModelMpi
    {
        private const int subdomainDataTag = 0;

        //TODO: This does not guarantee that the model also uses the same elementFactory for the elements of this process's 
        //      subdomain.
        private readonly IXFiniteElementFactory elementFactory;

        public XModelMpiCentralized(ProcessDistribution processDistribution, Func<XModel> createModel,
            IXFiniteElementFactory elementFactory) : base(processDistribution)
        {
            this.elementFactory = elementFactory;
            if (processDistribution.IsMasterProcess) this.model = createModel();
            else
            {
                this.model = new XModel();
                foreach (int s in processDistribution.GetSubdomainIDsOfProcess(processDistribution.OwnRank))
                {
                    this.model.Subdomains[s] = new XSubdomain(s);
                }
            }
        }

        public IDomain2DBoundary Boundary => this.model.Boundary;

        public new IDofSerializer DofSerializer
        {
            set => model.DofSerializer = value;
        }

        public XModel RawModel => model;

        public XSubdomain GetXSubdomain(int subdomainID)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomainID);
            return model.Subdomains[subdomainID];
        }

        protected override void ScatterSubdomainData()
        {
            HashSet<int> allSubdomainIDs = null;
            if (procs.IsMasterProcess) allSubdomainIDs = new HashSet<int>(model.EnumerateSubdomains().Select(sub => sub.ID));
            ScatterSubdomains(allSubdomainIDs);
        }
        
        public void ScatterSubdomains(HashSet<int> modifiedSubdomains)
        {
            // Broadcast which subdomains are modified
            int[] subdomainIDs = null;
            if (procs.IsMasterProcess) subdomainIDs = modifiedSubdomains.ToArray();
            MpiUtilities.BroadcastArray(procs.Communicator, ref subdomainIDs, procs.MasterProcess);
            if (!procs.IsMasterProcess) modifiedSubdomains = new HashSet<int>(subdomainIDs);

            // Scatter the modified subdomain data
            var activeSubdomains = new ActiveSubdomains(procs, s => modifiedSubdomains.Contains(s));
            var transferrer = new TransferrerPerSubdomain(procs);
            PackSubdomainData<XSubdomain, XSubdomainDto> packData =
                (s, subdomain) => XSubdomainDto.Serialize(subdomain, model.DofSerializer);
            UnpackSubdomainData<XSubdomain, XSubdomainDto> unpackData = 
                (s, subdomainDto) => subdomainDto.Deserialize(model.DofSerializer, elementFactory);
            Dictionary<int, XSubdomain> subdomainsOfProcess = transferrer.ScatterToSomeSubdomainsPacked(
                model.Subdomains, packData, unpackData, activeSubdomains);

            // Connect data structures of each subdomain in slave processes
            if (!procs.IsMasterProcess)
            {
                foreach (int s in subdomainsOfProcess.Keys)
                {
                    XSubdomain subdomain = model.Subdomains[s];
                    subdomain.ReplaceEntitiesWith(subdomainsOfProcess[s]);
                    subdomain.ConnectDataStructures();
                }
            }
        }

        public void ScatterSubdomainsState()
        {
            ScatterSubdomainsState(sub => sub.ConnectivityModified, (sub, modified) => sub.ConnectivityModified = modified);
            ScatterSubdomainsState(sub => sub.StiffnessModified, (sub, modified) => sub.StiffnessModified = modified);
        }

        private void ScatterSubdomainsState(Func<ISubdomain, bool> inquireStateModified, 
            Action<ISubdomain, bool> setStateModified)
        {
            // Prepare data in master
            Dictionary<int, bool> allSubdomainsState_master = null;
            if (procs.IsMasterProcess)
            {
                allSubdomainsState_master = new Dictionary<int, bool>();
                foreach (XSubdomain subdomain in model.Subdomains.Values)
                {
                    allSubdomainsState_master[subdomain.ID] = inquireStateModified(subdomain);
                }
            }

            // Scatter them to all processes
            var transferrer = new TransferrerAltogetherFlattened(procs);
            Dictionary<int, bool> processSubdomainsState = transferrer.ScatterToAllSubdomains(allSubdomainsState_master);

            // Update subdomains in other processes
            if (!procs.IsMasterProcess)
            {
                foreach (int s in processSubdomainsState.Keys)
                {
                    XSubdomain subdomain = model.Subdomains[s];
                    setStateModified(subdomain, processSubdomainsState[s]);
                }
            }
        }
    }
}
