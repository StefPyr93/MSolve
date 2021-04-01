using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using MPI;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer
{
    public static class TransferrerTestUtilities
    {
        private const int numProcesses = 4;

        public enum SubdomainDistribution
        {
            OnePerProcess, Uniform, Variable
        }

        public enum TransferrerChoice
        {
            PerSubdomain, AltogetherFlattened
        }

        public static ActiveSubdomains DetermineActiveSubdomains(ProcessDistribution procs)
        {
            // Every third subdomain is active
            return new ActiveSubdomains(procs, s => (s + 1) % 3 == 0);
        }

        public static ProcessDistribution DetermineProcesses(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            if (numProcesses != 4) throw new MpiProcessesException("Number of MPI processes must be exactly 4");
            int master = 0;
            var processesToSubdomains = new int[numProcesses][];
            int numSubdomains = 0;
            for (int p = 0; p < numProcesses; ++p)
            {
                int numSubdomainsOfThisProcess;
                if (subdomainDistribution == SubdomainDistribution.OnePerProcess) numSubdomainsOfThisProcess = 1;
                else if (subdomainDistribution == SubdomainDistribution.Uniform) numSubdomainsOfThisProcess = 5;
                else if (subdomainDistribution == SubdomainDistribution.Variable) numSubdomainsOfThisProcess = p + 1;
                else throw new NotImplementedException();
                processesToSubdomains[p] = new int[numSubdomainsOfThisProcess];
                {
                    for (int i = 0; i < numSubdomainsOfThisProcess; ++i)
                    {
                        processesToSubdomains[p][i] = numSubdomains++;
                    }
                }
            }
            return new ProcessDistribution(Communicator.world, master, processesToSubdomains);
        }

        public static ISubdomainDataTransferrer DetermineTransferrer(TransferrerChoice transferrerChoice, ProcessDistribution procs)
        {
            if (transferrerChoice == TransferrerChoice.PerSubdomain) return new TransferrerPerSubdomain(procs);
            else if (transferrerChoice == TransferrerChoice.AltogetherFlattened) return new TransferrerAltogetherFlattened(procs);
            else throw new NotImplementedException();
        }
    }
}
