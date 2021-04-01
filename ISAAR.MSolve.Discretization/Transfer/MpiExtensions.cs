using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;

namespace ISAAR.MSolve.Discretization.Transfer
{
    public static class MpiExtensions
    {
        public static Dictionary<ISubdomain, T> ChangeKey<T>(this Dictionary<int, T> subdomainIDsToData, IModel model)
        {
            var subdomainsToData = new Dictionary<ISubdomain, T>();
            foreach (var pair in subdomainIDsToData)
            {
                ISubdomain subdomain = model.GetSubdomain(pair.Key);
                subdomainsToData[subdomain] = pair.Value;
            }
            return subdomainsToData;
        }

        public static Dictionary<int, T> ChangeKey<T>(this Dictionary<ISubdomain, T> subdomainsToData)
        {
            var subdomainIDsToData = new Dictionary<int, T>();
            foreach (var pair in subdomainsToData) subdomainIDsToData[pair.Key.ID] = pair.Value;
            return subdomainIDsToData;
        }

        //TODO: Shouldn't this be accessed by model itself? But then what happens for master process, which contains all subdomains?
        public static Dictionary<int, ISubdomain> GetSubdomainsOfProcess(this ProcessDistribution procs, IModel model)
        {
            var processSubdomains = new Dictionary<int, ISubdomain>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank)) processSubdomains[s] = model.GetSubdomain(s);
            return processSubdomains;
        }
    }
}
