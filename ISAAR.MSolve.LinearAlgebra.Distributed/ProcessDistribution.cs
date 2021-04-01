using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed.Exceptions;
using MPI;

//TODO: If there are point2point communications between two processes, should this be transfered to all of them?
namespace ISAAR.MSolve.LinearAlgebra.Distributed
{
    /// <summary>
    /// This should be present in all processes.
    /// </summary>
    public class ProcessDistribution
    {
        private readonly int[][] processesToSubdomains;
        //private readonly Dictionary<int, int> subdomainsToProcesses;

        public ProcessDistribution(Intracommunicator comm, int masterProcess, int[][] processRanksToSubdomainIDs)
        {
            this.Communicator = comm;
            this.IsMasterProcess = comm.Rank == masterProcess;
            this.MasterProcess = masterProcess;
            this.OwnRank = comm.Rank;

            this.processesToSubdomains = processRanksToSubdomainIDs;
            //this.subdomainsToProcesses = new Dictionary<int, int>();
            NumSubdomainsTotal = 0;
            for (int p = 0; p < comm.Size; ++p)
            {
                NumSubdomainsTotal += processRanksToSubdomainIDs[p].Length;
                //foreach (int s in processRanksToSubdomainIDs[p]) this.subdomainsToProcesses[s] = p;
            }
        }

        public Intracommunicator Communicator { get; }
        public bool IsMasterProcess { get; }
        public int MasterProcess { get; }
        public int NumSubdomainsTotal { get; }
        public int OwnRank { get; }

        /// <summary>
        /// Example: 3 processes & 8 subdomains:
        /// process 0: s0, s1, s2
        /// process 1: s3, s4, s5
        /// process 2: s6, s7
        /// </summary>
        /// <param name="numProcesses"></param>
        /// <param name="numSubdomains"></param>
        public static ProcessDistribution CreateDistribution(int numProcesses, int numSubdomains)
        {
            if ((numProcesses < 2) || (numProcesses > numSubdomains))
            {
                throw new MpiProcessesException($"Number of MPI processes must belong to [2, {numSubdomains}]");
            }

            // Gather a set of subdomains that are multiple of the number of processes and distribute them evenly. 
            // Then each of the remainder subdomains goes to one of the first processes. 
            int div = numSubdomains / numProcesses;
            int mod = numSubdomains % numProcesses;
            var numSubdomainsPerProcess = new int[numProcesses];
            for (int p = 0; p < numProcesses; ++p) numSubdomainsPerProcess[p] = div;
            for (int m = 0; m < mod; ++m) ++numSubdomainsPerProcess[m];


            int subdomainID = 0; //TODO: Risky
            int[][] processesToSubdomains = new int[numProcesses][];
            for (int p = 0; p < numProcesses; ++p)
            {
                processesToSubdomains[p] = new int[numSubdomainsPerProcess[p]];
                for (int i = 0; i < numSubdomainsPerProcess[p]; ++i) processesToSubdomains[p][i] = subdomainID++;
            }

            int master = 0;
            var procs = new ProcessDistribution(MPI.Communicator.world, master, processesToSubdomains);
            //Console.WriteLine($"(process {procs.OwnRank}) Hello World!"); // Run this to check if MPI works correctly.
            //PrintProcessDistribution(procs);
            return procs;
        }

        [Conditional("DEBUG")]
        public void CheckProcessIsMaster()
        {
            if (!IsMasterProcess) throw new MpiException(
                $"Process {OwnRank}: Only defined for master process (rank = {MasterProcess})");
        }

        [Conditional("DEBUG")]
        public void CheckProcessMatchesSubdomain(int subdomainID)
        {
            bool isStored = processesToSubdomains[OwnRank].Contains(subdomainID);
            if (!isStored) throw new MpiException(
                $"Process {OwnRank}: This process does not have access to subdomain {subdomainID}");
        }

        [Conditional("DEBUG")]
        public void CheckProcessMatchesSubdomainUnlessMaster(int subdomainID)
        {
            if (IsMasterProcess) return;
            bool isStored = processesToSubdomains[OwnRank].Contains(subdomainID);
            if (!isStored) throw new MpiException(
                $"Process {OwnRank}: This process does not have access to subdomain {subdomainID}");
        }

        public int[] GetNumSubdomainsPerProcess()
        {
            var counts = new int[Communicator.Size];
            for (int p = 0; p < counts.Length; ++p) counts[p] = processesToSubdomains[p].Length;
            return counts;
        }

        /// <summary>
        /// The subdomain IDs are in the same order for all processes.
        /// </summary>
        /// <param name="processRank"></param>
        public int[] GetSubdomainIDsOfProcess(int processRank) => processesToSubdomains[processRank];

        public void PrintProcessDistribution()
        {
            if (IsMasterProcess)
            {
                Console.WriteLine("---------------------------------------");
                Console.WriteLine($"Printing from process {OwnRank}");
                for (int p = 0; p < Communicator.Size; ++p)
                {
                    var builder1 = new StringBuilder();
                    builder1.Append($"Process {p} - subdomain IDs: ");
                    foreach (int s in GetSubdomainIDsOfProcess(p))
                    {
                        builder1.Append(s + " ");
                    }
                    Console.WriteLine(builder1);
                }
            }

            if (IsMasterProcess)
            {
                Console.WriteLine("---------------------------------------");
                Console.WriteLine($"Printing from each process:");
            }
            Communicator.Barrier();
            var builder2 = new StringBuilder();
            builder2.Append($"Process {OwnRank} - subdomain IDs: ");
            foreach (int s in GetSubdomainIDsOfProcess(OwnRank))
            {
                builder2.Append(s + " ");
            }
            Console.WriteLine(builder2);
        }
    }
}
