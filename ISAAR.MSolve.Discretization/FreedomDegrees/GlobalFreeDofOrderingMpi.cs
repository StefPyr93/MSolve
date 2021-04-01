using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using MPI;

//TODO: For now each global-subdomain mappings are stored in both in master process and in the corresponding subdomain. However 
//      it may be more efficient to keep them in master and do the work there
//TODO: I implemented the methods as being executed simultaneously for all processes/subdomains. 
//      This is in contrast to the serial implementation that only does one subdomain at a time, meaning different schematics.
//      To use the interface polymorphically, the interface should specify both an all-together version and one-at-a-time.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class GlobalFreeDofOrderingMpi : IGlobalFreeDofOrdering
    {
        private readonly DofTable globalFreeDofs_master;
        private readonly IModel model;
        private readonly int numGlobalFreeDofs_master;
        private readonly ProcessDistribution procs;

        private bool hasGatheredSubdomainOrderings = false;
        private bool hasCreatedSubdomainGlobalMaps = false;
        private bool hasScatteredSubdomainGlobalMaps = false;

        /// <summary>
        /// Master contains all orderings. All other process only contain the corresponding subdomain data.
        /// </summary>
        private Dictionary<int, ISubdomainFreeDofOrdering> subdomainDofOrderings;

        /// <summary>
        /// Master contains all maps. All other process only contain the corresponding subdomain data.
        /// </summary>
        private Dictionary<int, int[]> subdomainToGlobalDofMaps;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="numGlobalFreeDofs">Will be ignored for every process other than <paramref name="processes.MasterProcess"/>.</param>
        /// <param name="globalFreeDofs">Will be ignored for every process other than <paramref name="processes.MasterProcess"/>.</param>
        public GlobalFreeDofOrderingMpi(ProcessDistribution processDistribution, IModel model, int numGlobalFreeDofs,
            DofTable globalFreeDofs)
        {
            this.procs = processDistribution;
            this.model = model;
            this.numGlobalFreeDofs_master = numGlobalFreeDofs;
            this.globalFreeDofs_master = globalFreeDofs;

            this.subdomainDofOrderings = new Dictionary<int, ISubdomainFreeDofOrdering>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomainDofOrderings[subdomain.ID] = subdomain.FreeDofOrdering;
            }
        }

        public DofTable GlobalFreeDofs
        {
            get
            {
                procs.CheckProcessIsMaster();
                return globalFreeDofs_master;
            }
        }

        public int NumGlobalFreeDofs //TODO: This can be broadcasted tbh
        {
            get
            {
                procs.CheckProcessIsMaster();
                return numGlobalFreeDofs_master;
            }
        }

        /// <summary>
        /// This method must only be called by the master process, after all data are in the correct memory spaces and 
        /// <see cref="CreateSubdomainGlobalMaps"/> has been called.
        /// </summary>
        /// <param name="subdomain"></param>
        /// <param name="subdomainVector">Each process has its own.</param>
        /// <param name="globalVector">Only exists in master process.</param>
        public void AddVectorSubdomainToGlobal(ISubdomain subdomain, IVectorView subdomainVector, IVector globalVector)
        {
            procs.CheckProcessIsMaster();
            int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain.ID];
            globalVector.AddIntoThisNonContiguouslyFrom(subdomainToGlobalDofs, subdomainVector);
        }

        public void AddVectorSubdomainToGlobalMeanValue(ISubdomain subdomain, IVectorView subdomainVector,
            IVector globalVector) => throw new NotImplementedException();

        public void CreateSubdomainGlobalMaps()
        {
            if (hasCreatedSubdomainGlobalMaps) return;

            GatherSubdomainDofOrderings();

            // Create and store all mapping arrays in master
            if (procs.IsMasterProcess)
            {
                subdomainToGlobalDofMaps =
                    GlobalFreeDofOrderingUtilities.CalcSubdomainGlobalMappings(globalFreeDofs_master, subdomainDofOrderings);
            }
            hasCreatedSubdomainGlobalMaps = true;
        }

        //TODO: this does not work if called only be master or for a subdomain that does not correspond to the process (even for master).
        /// <summary>
        /// 
        /// </summary>
        /// <param name="subdomain"></param>
        /// <param name="globalVector">Only exists in master process.</param>
        /// <param name="subdomainVector">Each process has its own.</param>
        public void ExtractVectorSubdomainFromGlobal(ISubdomain subdomain, IVectorView globalVector, IVector subdomainVector)
        {
            throw new NotImplementedException();
            
            // These were for the case when each process had only 1 subdomain.
            //procs.CheckProcessMatchesSubdomain(subdomain.ID);

            //ScatterSubdomainGlobalMaps();

            //// Broadcast globalVector
            ////TODO: The next is stupid, since it copies the vector to an array, while I could access its backing storage in 
            ////      most cases. I need a class that handles transfering the concrete vector class. That would live in an 
            ////      LinearAlgebra.MPI project
            //double[] globalArray = null;
            //if (procs.IsMasterProcess) globalArray = globalVector.CopyToArray();
            //MpiUtilities.BroadcastArray<double>(procs.Communicator, ref globalArray, procs.MasterProcess);
            //if (!procs.IsMasterProcess) globalVector = Vector.CreateFromArray(globalArray);

            //// Then do the actual work
            //int[] subdomainToGlobalDofs = subdomainToGlobalDofMaps[subdomain.ID];
            //subdomainVector.CopyNonContiguouslyFrom(globalVector, subdomainToGlobalDofs);
        }

        public void GatherSubdomainDofOrderings()
        {
            if (hasGatheredSubdomainOrderings) return;

            // Prepare data in each process
            var processOrderings = new Dictionary<int, DofTable>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                processOrderings[s] = subdomainDofOrderings[s].FreeDofs;
            }

            // Gather data in master
            var transferrer = new TransferrerPerSubdomain(procs);
            var tableSerializer = new DofTableSerializer(model.DofSerializer);
            GetArrayLengthOfPackedData<DofTable> getPackedDataLength = (s, table) => tableSerializer.CalcPackedLength(table);
            PackSubdomainDataIntoArray<DofTable, int> packData =
                (s, table, buffer, offset) => tableSerializer.PackTableIntoArray(table, buffer, offset);
            UnpackSubdomainDataFromArray<DofTable, int> unpackData =
                (s, buffer, start, end) => tableSerializer.UnpackTableFromArray(buffer, start, end, model.GetNode);
            Dictionary<int, DofTable> allOrderings = transferrer.GatherFromAllSubdomainsPacked(processOrderings,
                getPackedDataLength, packData, unpackData);

            // Store the received data in master
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.OwnRank) continue;
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        DofTable ordering = allOrderings[sub];
                        subdomainDofOrderings[sub] = new SubdomainFreeDofOrderingCaching(ordering.EntryCount, ordering);
                    }
                }
            }

            hasGatheredSubdomainOrderings = true;
            procs.Communicator.Barrier();
        }

        public ISubdomainFreeDofOrdering GetSubdomainDofOrdering(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomain.ID);
            return subdomainDofOrderings[subdomain.ID];
        }

        public int[] MapSubdomainToGlobalDofs(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomain.ID);
            ScatterSubdomainGlobalMaps();
            return subdomainToGlobalDofMaps[subdomain.ID];
        }

        private void ScatterSubdomainGlobalMaps()
        {
            if (hasScatteredSubdomainGlobalMaps) return;

            CreateSubdomainGlobalMaps();
            var transferrer = new TransferrerPerSubdomain(procs);
            Dictionary<int, int[]> processMappings = transferrer.ScatterToAllSubdomains(subdomainToGlobalDofMaps);

            // Store each one in processes other than master, since it already has all of them.
            if (!procs.IsMasterProcess) subdomainToGlobalDofMaps = processMappings;
            
            hasScatteredSubdomainGlobalMaps = true;
        }
    }
}
