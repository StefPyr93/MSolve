using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Transfer
{
    public class TransferrerPerSubdomain : ISubdomainDataTransferrer
    {
        private readonly ProcessDistribution procs;

        public TransferrerPerSubdomain(ProcessDistribution processDistribution)
        {
            this.procs = processDistribution;
        }

        public Dictionary<int, T> GatherFromAllSubdomains<T>(Dictionary<int, T> processSubdomainsData)
            => GatherFromAllSubdomainsPacked<T, T>(processSubdomainsData, (s, data) => data, (s, data) => data);

        public Dictionary<int, T[]> GatherFromAllSubdomains<T>(Dictionary<int, T[]> processSubdomainsData)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return GatherFromSomeSubdomains<T>(processSubdomainsData, activeSubdomains);
        }

        public Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return GatherFromSomeSubdomainsPacked<TRaw, TPacked>(processSubdomainsData, packData, unpackData, activeSubdomains);
        }

        public Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return GatherFromSomeSubdomainsPacked<TRaw, TPacked>(processSubdomainsData, getPackedDataLength, packData, 
                unpackData, activeSubdomains);
        }

        public Dictionary<int, T> GatherFromSomeSubdomains<T>(Dictionary<int, T> processSubdomainsData, 
            ActiveSubdomains activeSubdomains)
        {
            return GatherFromSomeSubdomainsPacked<T, T>(processSubdomainsData, (s, data) => data, (s, data) => data,
                activeSubdomains);
        }

        public Dictionary<int, T[]> GatherFromSomeSubdomains<T>(Dictionary<int, T[]> processSubdomainsData,
            ActiveSubdomains activeSubdomains)
        {
            // Pack and send the subdomain data to the corresponding process, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                var allSubdomainsData = new Dictionary<int, T[]>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) // Just copy the references
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s)) allSubdomainsData[s] = processSubdomainsData[s];
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                allSubdomainsData[s] = MpiUtilities.ReceiveArray<T>(procs.Communicator, p, s);
                            }
                        }
                    }
                }
                return allSubdomainsData;
            }
            else
            {
                // Send the data of each subdomain to master
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        MpiUtilities.SendArray<T>(procs.Communicator, processSubdomainsData[s], procs.MasterProcess, s);
                    }
                }
                return null;
            }
        }

        public Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData,
            ActiveSubdomains activeSubdomains)
        {
            // Gather and unpack the subdomain data from the corresponding processes, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                var allSubdomainsData = new Dictionary<int, TRaw>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) // Just copy the references
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s)) allSubdomainsData[s] = processSubdomainsData[s]; 
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                TPacked packed = procs.Communicator.Receive<TPacked>(p, s);
                                allSubdomainsData[s] = unpackData(s, packed);
                            }
                        }
                    }
                }
                return allSubdomainsData;
            }
            else
            {
                // At first pack all subdomain data of this process, so that master process doe not have to wait inbetween
                //TODO: For the first process, master must wait longer this way. Use asynchronous sends or another thread instead.
                int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                var processDataPacked = new Dictionary<int, TPacked>(processSubdomains.Length);
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processDataPacked[s] = packData(s, processSubdomainsData[s]);
                    }
                }

                // Then send the packed subdomain data to master
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        procs.Communicator.Send<TPacked>(processDataPacked[s], procs.MasterProcess, s);
                        processDataPacked.Remove(s); // Free up some memory
                    }
                }
                return null;
            }
        }


        public Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains)
        {
            Dictionary<int, TRaw> allSubdomainsData = null;

            // Pack and send the subdomain data to the corresponding process, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                allSubdomainsData = new Dictionary<int, TRaw>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) // Just copy the references
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s)) allSubdomainsData[s] = processSubdomainsData[s];
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                TPacked[] packed = MpiUtilities.ReceiveArray<TPacked>(procs.Communicator, p, s);
                                allSubdomainsData[s] = unpackData(s, packed, 0, packed.Length);
                            }
                        }
                    }
                }
                return allSubdomainsData;
            }
            else
            {
                // At first pack all subdomain data of this process, so that master process doe not have to wait inbetween
                //TODO: For the first process, master must wait longer this way. Use asynchronous sends or another thread instead.
                int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                var processDataPacked = new Dictionary<int, TPacked[]>(processSubdomains.Length);
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        TRaw rawData = processSubdomainsData[s];
                        int packedLength = getPackedDataLength(s, rawData);
                        var packedData = new TPacked[packedLength];
                        packData(s, rawData, packedData, 0);
                        processDataPacked[s] = packedData;
                    }
                }

                // Then send the packed subdomain data to master
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        MpiUtilities.SendArray<TPacked>(procs.Communicator, processDataPacked[s], procs.MasterProcess, s);
                        processDataPacked.Remove(s); // Free up some memory
                    }
                }
                return null;
            }
        }

        public Dictionary<int, T> ScatterToAllSubdomains<T>(Dictionary<int, T> allSubdomainsData_master)
            => ScatterToAllSubdomainsPacked<T, T>(allSubdomainsData_master, (s, data) => data, (s, data) => data);

        public Dictionary<int, T[]> ScatterToAllSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return ScatterToSomeSubdomains<T>(allSubdomainsData_master, activeSubdomains);
        }

        /// <summary>
        /// This method returns null in master process. For other processes, it returns a Dictionary with the data for each 
        /// associated subdomain.
        /// </summary>
        public Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return ScatterToSomeSubdomainsPacked<TRaw, TPacked>(allSubdomainsData_master, packData, unpackData, activeSubdomains);
        }

        public Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData)
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => true);
            return ScatterToSomeSubdomainsPacked<TRaw, TPacked>(allSubdomainsData_master, getPackedDataLength, packData, 
                unpackData, activeSubdomains);
        }

        /// <summary>
        /// This method returns null in master process. For other processes, it returns a Dictionary with the data for each 
        /// associated subdomain or an empty Dictionary if none of that process's subdomains are active.
        /// </summary>
        public Dictionary<int, T> ScatterToSomeSubdomains<T>(Dictionary<int, T> allSubdomainsData_master, 
            ActiveSubdomains activeSubdomains)
        {
            return ScatterToSomeSubdomainsPacked<T, T>(allSubdomainsData_master, (s, data) => data, (s, data) => data,
                activeSubdomains);
        }

        public Dictionary<int, T[]> ScatterToSomeSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master,
            ActiveSubdomains activeSubdomains)
        {
            int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            var processData = new Dictionary<int, T[]>(processSubdomains.Length);

            // Pack send the subdomain data to the corresponding process, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess)
                    {
                        foreach (int s in processSubdomains)
                        {
                            if (activeSubdomains.IsActive(s)) processData[s] = allSubdomainsData_master[s];
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                T[] data = allSubdomainsData_master[s];
                                MpiUtilities.SendArray<T>(procs.Communicator, data, p, s);
                            }
                        }
                    }
                }
            }
            else
            {
                // Receive all subdomain data of this process
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processData[s] = MpiUtilities.ReceiveArray<T>(procs.Communicator, procs.MasterProcess, s);
                    }
                }
            }

            return processData;
        }

        public Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData, 
            ActiveSubdomains activeSubdomains)
        {
            int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            var processData = new Dictionary<int, TRaw>(processSubdomains.Length);

            // Pack and send the subdomain data to the corresponding process, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess)
                    {
                        foreach (int s in processSubdomains)
                        {
                            if (activeSubdomains.IsActive(s)) processData[s] = allSubdomainsData_master[s];
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                TPacked packed = packData(s, allSubdomainsData_master[s]);
                                procs.Communicator.Send<TPacked>(packed, p, s);
                            }
                        }
                    }
                }
            }
            else
            {
                // At first, receive all subdomain data of this process, so that master can continue to the next process.
                var processDataPacked = new Dictionary<int, TPacked>(processSubdomains.Length);
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processDataPacked[s] = procs.Communicator.Receive<TPacked>(procs.MasterProcess, s);
                    }
                }

                // Then unpack and return the subdomain data in each process
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processData[s] = unpackData(s, processDataPacked[s]);
                        processDataPacked.Remove(s); // Free up some memory
                    }
                }
            }

            return processData;
        }

        public Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains)
        {
            int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            var processData = new Dictionary<int, TRaw>(processSubdomains.Length);

            // Pack and send the subdomain data to the corresponding process, one subdomain at a time
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess)
                    {
                        foreach (int s in processSubdomains)
                        {
                            if (activeSubdomains.IsActive(s)) processData[s] = allSubdomainsData_master[s];
                        }
                    }
                    else
                    {
                        foreach (int s in procs.GetSubdomainIDsOfProcess(p))
                        {
                            if (activeSubdomains.IsActive(s))
                            {
                                TRaw rawData = allSubdomainsData_master[s];
                                int packedLength = getPackedDataLength(s, rawData);
                                var packedData = new TPacked[packedLength];
                                packData(s, rawData, packedData, 0);
                                MpiUtilities.SendArray<TPacked>(procs.Communicator, packedData, p, s);
                            }
                        }
                    }
                }
            }
            else
            {
                // At first, receive all subdomain data of this process, so that master can continue to the next process.
                var processDataPacked = new Dictionary<int, TPacked[]>(processSubdomains.Length);
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processDataPacked[s] = MpiUtilities.ReceiveArray<TPacked>(procs.Communicator, procs.MasterProcess, s);
                    }
                }

                // Then unpack and return the subdomain data in each process
                foreach (int s in processSubdomains)
                {
                    if (activeSubdomains.IsActive(s))
                    {
                        processData[s] = unpackData(s, processDataPacked[s], 0, processDataPacked[s].Length);
                        processDataPacked.Remove(s); // Free up some memory
                    }
                }
            }

            return processData;
        }
    }
}
