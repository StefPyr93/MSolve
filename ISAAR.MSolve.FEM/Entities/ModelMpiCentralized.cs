using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Transfer;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;

namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// Master process stores all data. Other processes store only data of their respective subdomains. Master process creates 
    /// the whole model and then scatters the subdomain data to their respective processes. Unfortunately this scattering 
    /// currently does not work as expected. If a process contains more than one subdomains, their common boundary nodes are
    /// created twice, which then causes problems in the DDM solver.
    /// </summary>
    public class ModelMpiCentralized : ModelMpiCentralizedBase<Model>
    {
        public ModelMpiCentralized(ProcessDistribution processDistribution, Func<Model> createModel) : 
            base(processDistribution)
        {
            if (processDistribution.IsMasterProcess) this.model = createModel();
            else this.model = new Model();
        }

        //TODO: This does not work as expected. If a process contains more than one subdomains, their common boundary nodes are
        // created twice, which then causes problems in the DDM solver.
        protected override void ScatterSubdomainData()
        {
            // Scatter subdomain data to all processes
            var transferrer = new TransferrerPerSubdomain(procs);
            PackSubdomainData<Subdomain, SubdomainDto> packData = 
                (id, subdomain) => SubdomainDto.Serialize(subdomain, DofSerializer);
            UnpackSubdomainData<Subdomain, SubdomainDto> unpackData =
                (id, subdomainDto) => subdomainDto.Deserialize(DofSerializer);
            Dictionary<int, Subdomain> subdomainsOfProcess = transferrer.ScatterToAllSubdomainsPacked(
                model.SubdomainsDictionary, packData, unpackData);

            if (!procs.IsMasterProcess)
            {
                // Add the subdomains to the model
                foreach (Subdomain subdomain in subdomainsOfProcess.Values)
                {
                    model.SubdomainsDictionary[subdomain.ID] = subdomain;
                    subdomain.ConnectDataStructures();
                }
            }
        }
    }
}
