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
    public abstract class ModelMpiCentralizedBase<TModel> : IModelMpi
        where TModel : IModel
    {
        protected readonly Dictionary<int, Cluster> clusters = new Dictionary<int, Cluster>();
        protected readonly ProcessDistribution procs;
        protected TModel model; // This is set in the concrete class's contructor for master and after scattering the subdomains for the other processes.

        protected ModelMpiCentralizedBase(ProcessDistribution processDistribution)
        {
            this.procs = processDistribution;
        }

        public Dictionary<int, Cluster> Clusters
        {
            get
            {
                procs.CheckProcessIsMaster();
                return clusters;
            }
        }

        public Table<INode, IDofType, double> Constraints
        {
            get
            {
                procs.CheckProcessIsMaster();
                return model.Constraints;
            }
        }

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

        public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads
        {
            get
            {
                procs.CheckProcessIsMaster();
                return model.MassAccelerationHistoryLoads;
            }
        }

        public int NumClusters
        {
            get
            {
                procs.CheckProcessIsMaster();
                return clusters.Count;
            }
        }

        public int NumElements
        {
            get
            {
                procs.CheckProcessIsMaster();
                return model.NumElements;
            }
        }

        public int NumNodes
        {
            get
            {
                procs.CheckProcessIsMaster();
                return model.NumNodes;
            }
        }

        public int NumSubdomains
        {
            get
            {
                procs.CheckProcessIsMaster();
                return model.NumSubdomains;
            }
        }

        public void ApplyLoads()
        {
            //model.ApplyLoads(); //TODO: This does not work in MPI environment.
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomain.Forces.Clear();
            }
        }

        public void ApplyMassAccelerationHistoryLoads(int timeStep)
        {
            procs.CheckProcessIsMaster();
            model.ApplyMassAccelerationHistoryLoads(timeStep);
        }

        public void ConnectDataStructures()
        {
            if (procs.IsMasterProcess)
            {
                model.ConnectDataStructures();
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    var cluster = new Cluster(p);
                    foreach (int s in procs.GetSubdomainIDsOfProcess(p)) cluster.Subdomains.Add(model.GetSubdomain(s));
                    clusters[p] = cluster;
                }
            }
            // If it is not master, then just return
        }

        public IEnumerable<Cluster> EnumerateClusters()
        {
            return clusters.Values;
        }

        public IEnumerable<IElement> EnumerateElements()
        {
            procs.CheckProcessIsMaster();
            return model.EnumerateElements();
        }

        public IEnumerable<INode> EnumerateNodes()
        {
            procs.CheckProcessIsMaster();
            return model.EnumerateNodes();
        }

        public IEnumerable<ISubdomain> EnumerateSubdomains()
        {
            return model.EnumerateSubdomains();
        }

        public Cluster GetCluster(int clusterID)
        {
            //procs.CheckProcessMatchesClusterUnlessMaster(clusterID);
            return clusters[clusterID];
        }

        public IElement GetElement(int elementID)
        {
            procs.CheckProcessIsMaster();
            return model.GetElement(elementID);
        }

        public INode GetNode(int nodeID)
        {
            procs.CheckProcessIsMaster();
            return model.GetNode(nodeID);
        }

        public ISubdomain GetSubdomain(int subdomainID)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomainID);
            return model.GetSubdomain(subdomainID);
        }

        public void ScatterSubdomains()
        {
            ScatterSubdomainData();
            if (!procs.IsMasterProcess)
            {
                var cluster = new Cluster(procs.OwnRank);
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank)) cluster.Subdomains.Add(model.GetSubdomain(s));
                clusters[procs.OwnRank] = cluster;
            }
        }

        protected abstract void ScatterSubdomainData();
    }
}
