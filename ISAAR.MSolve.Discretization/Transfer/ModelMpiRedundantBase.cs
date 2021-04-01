using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;

//TODO: The redirection to IModel is necessary to safeguard against a process trying to access unavailable data, but it slows 
//      access to the items needed by client code. This is especially pronounced for GetNode(), GetElement(), etc that will be 
//      called a lot of times in succession. Can it be optimized?
//TODO: Can the boilerplate code be reduced without sacrificing performance?
//TODO: Isn't it wasteful to create all this indirection, when the only things actually implemented by this are IDofSerializer, 
//      ScatterSubdomains() and GetSubdomain()? Must I use polymorphism for them? On the other hand, once updating the model is
//      is considered, it will be nice to let ModelMpi handle the communication before and after updating model fields.
namespace ISAAR.MSolve.Discretization.Transfer
{
    public class ModelMpiRedundantBase<TModel> : IModelMpi
        where TModel : IModel
    {
        protected readonly ProcessDistribution procs;
        protected readonly TModel model;

        protected ModelMpiRedundantBase(ProcessDistribution processDistribution, Func<TModel> createModel)
        {
            this.procs = processDistribution;
            this.model = createModel(); // Create the whole model in all processes.
        }

        public Dictionary<int, Cluster> Clusters { get; } = new Dictionary<int, Cluster>();

        public Table<INode, IDofType, double> Constraints => model.Constraints;

        public IDofSerializer DofSerializer
        {
            get => model.DofSerializer;
            set => model.DofSerializer = value;
        }

        public IGlobalFreeDofOrdering GlobalDofOrdering
        {
            get => model.GlobalDofOrdering;
            set => model.GlobalDofOrdering = value;
        }

        public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads => model.MassAccelerationHistoryLoads;
        public int NumClusters => Clusters.Count;

        public int NumElements => model.NumElements;

        public int NumNodes => model.NumNodes;

        public int NumSubdomains => model.NumSubdomains;

        public void ApplyLoads()
        {
            //model.ApplyLoads(); //TODO: This does not work in MPI environment.
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomain.Forces.Clear();
            }
        }

        public void ApplyMassAccelerationHistoryLoads(int timeStep) => model.ApplyMassAccelerationHistoryLoads(timeStep);

        public void ConnectDataStructures()
        {
            model.ConnectDataStructures();
            for (int p = 0; p < procs.Communicator.Size; ++p)
            {
                var cluster = new Cluster(p);
                foreach (int s in procs.GetSubdomainIDsOfProcess(p)) cluster.Subdomains.Add(model.GetSubdomain(s));
                Clusters[p] = cluster;
            }
        }

        public IEnumerable<Cluster> EnumerateClusters() => Clusters.Values;

        public IEnumerable<IElement> EnumerateElements() => model.EnumerateElements();

        public IEnumerable<INode> EnumerateNodes() => model.EnumerateNodes();

        public IEnumerable<ISubdomain> EnumerateSubdomains() => model.EnumerateSubdomains();

        public Cluster GetCluster(int clusterID) => Clusters[clusterID];

        public IElement GetElement(int elementID) => model.GetElement(elementID);

        public INode GetNode(int nodeID) => model.GetNode(nodeID);

        public ISubdomain GetSubdomain(int subdomainID) => model.GetSubdomain(subdomainID);

        public void ScatterSubdomains() { } // Do nothing
    }
}
