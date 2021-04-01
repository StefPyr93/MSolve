using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.Transfer;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IModelMpi : IModel
    {
        //TODO: DofSerializer belongs here instead of IModel. However that would create a lot of trouble with strategies that 
        //      expect IModel due to implementing general interfaces and have MPI implementations that need IDofSerializer
        //IDofSerializer DofSerializer { get; }

        int NumClusters { get; }

        IEnumerable<Cluster> EnumerateClusters();

        Cluster GetCluster(int clusterID);

        void ScatterSubdomains(); //TODO: perhaps this should be incorporated in ConnectDataStructures()
    }
}
