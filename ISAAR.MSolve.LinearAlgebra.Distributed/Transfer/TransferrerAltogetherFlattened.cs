using System;
using System.Collections.Generic;
using System.Text;

//TODO: This class uses the flattened versions of MPI calls. Add another one that uses the regular ones and benchmark.
//TODO: I should avoid packing/unpacking the subdomains of the master process. Also avoid placing them onto arrays.
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Transfer
{
    public class TransferrerAltogetherFlattened : ISubdomainDataTransferrer
    {
        private readonly ProcessDistribution procs;

        public TransferrerAltogetherFlattened(ProcessDistribution processDistribution)
        {
            this.procs = processDistribution;
        }

        public Dictionary<int, T> GatherFromAllSubdomains<T>(Dictionary<int, T> processSubdomainsData)
            => GatherFromAllSubdomainsPacked<T, T>(processSubdomainsData, (s, data) => data, (s, data) => data);

        public Dictionary<int, T[]> GatherFromAllSubdomains<T>(Dictionary<int, T[]> processSubdomainsData)
        {
            GetArrayLengthOfPackedData<T[]> getPackedDataLength = (s, subdomainArray) => subdomainArray.Length;
            PackSubdomainDataIntoArray<T[], T> packData = (s, subdomainArray, buffer, offset)
                => Array.Copy(subdomainArray, 0, buffer, offset, subdomainArray.Length);
            UnpackSubdomainDataFromArray<T[], T> unpackData = (s, buffer, start, end) =>
            {
                var subdomainArray = new T[end - start];
                Array.Copy(buffer, start, subdomainArray, 0, end - start);
                return subdomainArray;
            };
            return GatherFromAllSubdomainsPacked(processSubdomainsData, getPackedDataLength, packData, unpackData);
        }

        public Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData)
        {
            // Pack the subdomain data in each process
            int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            var processDataPacked = new TPacked[processSubdomains.Length];
            for (int s = 0; s < processSubdomains.Length; ++s)
            {
                int sub = processSubdomains[s];
                TPacked packed = packData(sub, processSubdomainsData[sub]);
                processDataPacked[s] = packed;
            }

            // Gather all subdomain data in master
            int[] numSubdomainsPerProcess = procs.GetNumSubdomainsPerProcess();
            TPacked[] allDataPacked = procs.Communicator.GatherFlattened<TPacked>(
                processDataPacked, numSubdomainsPerProcess, procs.MasterProcess);

            // Unpack all subdomain data in master
            Dictionary<int, TRaw> allData = null;
            int offset = 0;
            if (procs.IsMasterProcess)
            {
                allData = new Dictionary<int, TRaw>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        TPacked packed = allDataPacked[offset];
                        allData[sub] = unpackData(sub, packed);
                        allDataPacked[offset] = default; // Free up some memory
                        ++offset;
                    }
                }
                
            }
            return allData;
        }

        public Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData)
        {
            // Determine the size of each subdomain's array in each process
            var processArraySizes = new Dictionary<int, int>();
            int processArraySizeTotal = 0;
            int[] processSubdomains = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            foreach (int sub in processSubdomains)
            {
                TRaw rawData = processSubdomainsData[sub];
                int packedSize = getPackedDataLength(sub, rawData);
                processArraySizes[sub] = packedSize;
                processArraySizeTotal += packedSize;
            }

            // Determine the size of each subdomain's packed data in master
            Dictionary<int, int> allPackedSizes = GatherFromAllSubdomains<int>(processArraySizes);
            int totalPackedSize = 0;
            int[] processTotalPackedSizes = null; 
            if (procs.IsMasterProcess)
            {
                processTotalPackedSizes = new int[procs.Communicator.Size];
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        int size = allPackedSizes[sub];
                        processTotalPackedSizes[p] += size;
                        totalPackedSize += size;
                    }
                }
            }

            // Pack the subdomain data in each process
            var processDataPacked = new TPacked[processArraySizeTotal];
            int offset = 0;
            for (int s = 0; s < processSubdomains.Length; ++s)
            {
                int sub = processSubdomains[s];
                packData(sub, processSubdomainsData[sub], processDataPacked, offset);
                offset += processArraySizes[sub];
            }

            // Gather all subdomain data in master
            TPacked[] allDataPacked = procs.Communicator.GatherFlattened<TPacked>(
                processDataPacked, processTotalPackedSizes, procs.MasterProcess);

            // Unpack all subdomain data in master
            Dictionary<int, TRaw> allData = null;
            int start = 0;
            if (procs.IsMasterProcess)
            {
                allData = new Dictionary<int, TRaw>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        int end = start + allPackedSizes[sub];
                        allData[sub] = unpackData(sub, allDataPacked, start, end);
                        start = end;
                    }
                }
            }
            return allData;
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
            throw new NotImplementedException();
        }

        public Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData, 
            ActiveSubdomains activeSubdomains)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, T> ScatterToAllSubdomains<T>(Dictionary<int, T> allSubdomainsData_master)
            => ScatterToAllSubdomainsPacked<T, T>(allSubdomainsData_master, (s, data) => data, (s, data) => data);

        public Dictionary<int, T[]> ScatterToAllSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master)
        {
            GetArrayLengthOfPackedData<T[]> getPackedDataLength = (s, subdomainArray) => subdomainArray.Length;
            PackSubdomainDataIntoArray<T[], T> packData = (s, subdomainArray, buffer, offset) 
                => Array.Copy(subdomainArray, 0, buffer, offset, subdomainArray.Length);
            UnpackSubdomainDataFromArray<T[], T> unpackData = (s, buffer, start, end) =>
            {
                var subdomainArray = new T[end - start];
                Array.Copy(buffer, start, subdomainArray, 0, end - start);
                return subdomainArray;
            };
            return ScatterToAllSubdomainsPacked(allSubdomainsData_master, getPackedDataLength, packData, unpackData);
        }

        public Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master, 
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData)
        {
            // Put all data in a global array
            TPacked[] allDataPacked = null;
            if (procs.IsMasterProcess)
            {
                allDataPacked = new TPacked[procs.NumSubdomainsTotal];
                int offset = 0;
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        TPacked packed = packData(sub, allSubdomainsData_master[sub]);
                        allDataPacked[offset++] = packed;
                    }
                }
            }

            // Scatter the global array
            int[] numSubdomainsPerProcess = procs.GetNumSubdomainsPerProcess();
            TPacked[] processDataPacked = procs.Communicator.ScatterFromFlattened(
                allDataPacked, numSubdomainsPerProcess, procs.MasterProcess);

            // Unpack the subdomain data of this process. In master only copy the references
            int[] subdomainsOfProcess = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            Dictionary<int, TRaw> processData = new Dictionary<int, TRaw>(subdomainsOfProcess.Length);
            if (procs.IsMasterProcess)
            {
                foreach (int sub in subdomainsOfProcess)
                {
                    processData[sub] = allSubdomainsData_master[sub];
                }
            }
            {
                for (int s = 0; s < subdomainsOfProcess.Length; ++s)
                {
                    int sub = subdomainsOfProcess[s];
                    processData[sub] = unpackData(sub, processDataPacked[s]);
                }
            }
            return processData;
        }

        public Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData)
        {
            // Determine the size of each subdomain's packed data in master
            Dictionary<int, int> allPackedSizes = null;
            int totalPackedSize = 0;
            var processTotalPackedSizes = new int[procs.Communicator.Size];
            if (procs.IsMasterProcess)
            {
                allPackedSizes = new Dictionary<int, int>();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        TRaw rawData = allSubdomainsData_master[sub];
                        int packedSize = getPackedDataLength(sub, rawData);
                        totalPackedSize += packedSize;
                        processTotalPackedSizes[p] += packedSize;
                        allPackedSizes[sub] = packedSize;
                    }
                }
            }

            // Determine the size of each subdomain's packed data in other processes
            Dictionary<int, int> processPackedSizes = ScatterToAllSubdomains<int>(allPackedSizes);
            if (!procs.IsMasterProcess)
            {
                foreach (int size in processPackedSizes.Values) processTotalPackedSizes[procs.OwnRank] += size;
            }

            // Put all data in a global array in master
            TPacked[] allDataPacked = null;
            if (procs.IsMasterProcess)
            {
                allDataPacked = new TPacked[totalPackedSize];
                int offset = 0;
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int sub in procs.GetSubdomainIDsOfProcess(p))
                    {
                        TRaw rawData = allSubdomainsData_master[sub];
                        packData(sub, rawData, allDataPacked, offset);
                        offset += getPackedDataLength(sub, rawData);
                    }
                }
            }

            // Scatter the global array
            // In processes other than master, the "counts" parameter is an int array that contains the actual count/size only 
            // at the index that corresponds to this process. All other entries are 0 and ignored.
            TPacked[] processDataPacked = procs.Communicator.ScatterFromFlattened(
                allDataPacked, processTotalPackedSizes, procs.MasterProcess);

            // Unpack the subdomain data of this process. In master only copy the references
            int[] subdomainsOfProcess = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
            Dictionary<int, TRaw> processData = new Dictionary<int, TRaw>(subdomainsOfProcess.Length);
            if (procs.IsMasterProcess)
            {
                foreach (int sub in subdomainsOfProcess)
                {
                    processData[sub] = allSubdomainsData_master[sub];
                }
            }
            else
            {
                int start = 0;
                for (int s = 0; s < subdomainsOfProcess.Length; ++s)
                {
                    int sub = subdomainsOfProcess[s];
                    int packedSize = processPackedSizes[sub];
                    int end = start + packedSize;
                    processData[sub] = unpackData(sub, processDataPacked, start, end);
                    start = end;
                }
            }
            return processData;
        }

        public Dictionary<int, T> ScatterToSomeSubdomains<T>(Dictionary<int, T> allSubdomainsData_master, 
            ActiveSubdomains activeSubdomains)
        {
            return ScatterToSomeSubdomainsPacked<T, T>(allSubdomainsData_master, (s, data) => data, (s, data) => data,
                activeSubdomains);
        }

        public Dictionary<int, T[]> ScatterToSomeSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master, 
            ActiveSubdomains activeSubdomains)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData, 
            ActiveSubdomains activeSubdomains)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master, 
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains)
        {
            throw new NotImplementedException();
        }
    }
}
