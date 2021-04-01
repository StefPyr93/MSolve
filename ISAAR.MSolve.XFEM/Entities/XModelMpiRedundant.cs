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
    public class XModelMpiRedundant : ModelMpiRedundantBase<XModel>, IXModelMpi
    {
        public XModelMpiRedundant(ProcessDistribution processDistribution, Func<XModel> createModel) : 
            base(processDistribution, createModel)
        {
        }

        public IDomain2DBoundary Boundary => this.model.Boundary;

        public XModel RawModel => model;

        public XSubdomain GetXSubdomain(int subdomainID)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomainID);
            return model.Subdomains[subdomainID];
        }

        public void ScatterSubdomains(HashSet<int> modifiedSubdomains)
        {
            // Do nothing
        }

        public void ScatterSubdomainsState()
        {
            // Do nothing
            //ScatterSubdomainsState(sub => sub.ConnectivityModified, (sub, modified) => sub.ConnectivityModified = modified);
            //ScatterSubdomainsState(sub => sub.StiffnessModified, (sub, modified) => sub.StiffnessModified = modified);
        }

        //private void ScatterSubdomainsState(Func<ISubdomain, bool> inquireStateModified,
        //    Action<ISubdomain, bool> setStateModified)
        //{
        //    // Prepare data in master
        //    Dictionary<int, bool> allSubdomainsState_master = null;
        //    if (procs.IsMasterProcess)
        //    {
        //        allSubdomainsState_master = new Dictionary<int, bool>();
        //        foreach (XSubdomain subdomain in model.Subdomains.Values)
        //        {
        //            allSubdomainsState_master[subdomain.ID] = inquireStateModified(subdomain);
        //        }
        //    }

        //    // Scatter them to all processes
        //    var transferrer = new TransferrerAltogetherFlattened(procs);
        //    Dictionary<int, bool> processSubdomainsState = transferrer.ScatterToAllSubdomains(allSubdomainsState_master);

        //    // Update subdomains in other processes
        //    if (!procs.IsMasterProcess)
        //    {
        //        foreach (int s in processSubdomainsState.Keys)
        //        {
        //            XSubdomain subdomain = model.Subdomains[s];
        //            setStateModified(subdomain, processSubdomainsState[s]);
        //        }
        //    }
        //}
    }
}
