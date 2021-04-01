using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPSubdomainGlobalMappingMpi
    {
        private readonly IStiffnessDistribution distribution;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public FetiDPSubdomainGlobalMappingMpi(ProcessDistribution processDistribution, IModel model,
            IFetiDPDofSeparator dofSeparator, IStiffnessDistribution distribution)
        {
            this.procs = processDistribution;
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.distribution = distribution;
        }

        public double CalcGlobalForcesNorm(Func<ISubdomain, Vector> getSubdomainForces)
        {
            //TODO: This can be optimized: calculate the dot product f*f for the internal dofs of each subdomain separately,
            //      only assemble global vector for the boundary dofs, find its dot product with itself, add the contributions
            //      for the internal dofs and finally apply SQRT(). This would greatly reduce the communication requirements.
            //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
            //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

            Vector globalForces = AssembleSubdomainVectors(getSubdomainForces);
            double norm = double.NaN;
            if (procs.IsMasterProcess) norm = globalForces.Norm2();
            procs.Communicator.Broadcast(ref norm, procs.MasterProcess); //TODO: Not sure if this is needed.
            return norm;
        }

        public Vector GatherGlobalDisplacements(Func<ISubdomain, Vector> getSubdomainFreeDisplacements)
        {
            return AssembleSubdomainVectors(sub =>
            {
                Vector u = getSubdomainFreeDisplacements(sub);
                FetiDPSubdomainGlobalMappingUtilities.ScaleSubdomainFreeDisplacements(dofSeparator, distribution, sub, u);
                return u;
            });
        }

        private Vector AssembleSubdomainVectors(Func<ISubdomain, Vector> getSubdomainVector)
        {
            //TODO: Subdomain vectors should not be gathered in master altogether. Instead the vector of each CLUSTER should be 
            //      received and then added

            //TODO: This class and the following one need serious refactoring
            ((GlobalFreeDofOrderingMpi)(model.GlobalDofOrdering)).CreateSubdomainGlobalMaps();

            // Prepare the subdomain vectors of each process
            var processVectors = new Dictionary<int, Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank)) //ERROR: For master this means that only its own subdomains are processes.
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                processVectors[s] = getSubdomainVector(subdomain);
            }

            // Set up packing/unpacking of vectors
            ////TODO: The next is stupid, since it copies the vector to an array, while I could access its backing storage in 
            ////      most cases. I need a class that handles transfering the concrete vector class. That would live in an 
            ////      LinearAlgebra.MPI project
            GetArrayLengthOfPackedData<Vector> getPackedDataLength = (s, vector) => vector.Length;
            PackSubdomainDataIntoArray<Vector, double> packData = (s, vector, buffer, offset) 
                => Array.Copy(vector.CopyToArray(), 0, buffer, offset, vector.Length);
            UnpackSubdomainDataFromArray<Vector, double> unpackData = (s, buffer, start, end) =>
            {
                var subdomainArray = new double[end - start];
                Array.Copy(buffer, start, subdomainArray, 0, end - start);
                return Vector.CreateFromArray(subdomainArray);
            };

            // Gather subdomain vectors in master 
            var transferrer = new TransferrerPerSubdomain(procs);
            Dictionary<int, Vector> allVectors = transferrer.GatherFromAllSubdomainsPacked(processVectors,
                getPackedDataLength, packData, unpackData);

            // Assemble the subdomain vectors into a global one
            Vector globalVector = null;
            if (procs.IsMasterProcess)
            {
                globalVector = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    model.GlobalDofOrdering.AddVectorSubdomainToGlobal(subdomain, allVectors[subdomain.ID], globalVector);
                }
            }
            return globalVector;
        }
    }
}
