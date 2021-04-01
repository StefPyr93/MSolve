using System;
using System.Collections.Generic;
using System.Text;

//TODO: Gradually replace ISubdomain.StiffnessModified and ISubdomain.ConnectivityModified with this
//TODO: So far this is just a wrapped Dictionary
namespace ISAAR.MSolve.LinearAlgebra.Distributed
{
    public class ActiveSubdomains
    {
        private readonly Dictionary<int, bool> areActive;

        public ActiveSubdomains(ProcessDistribution procs, Func<int, bool> querySubdomainActivation)
        {
            areActive = new Dictionary<int, bool>();
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    foreach (int s in procs.GetSubdomainIDsOfProcess(p)) areActive[s] = querySubdomainActivation(s);
                }
            }
            else
            {
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank)) areActive[s] = querySubdomainActivation(s);
            }
        }

        public bool IsActive(int subdomainID) => areActive[subdomainID];
    }
}
