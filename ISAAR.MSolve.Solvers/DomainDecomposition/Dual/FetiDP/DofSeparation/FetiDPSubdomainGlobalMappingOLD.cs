using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPSubdomainGlobalMappingOLD
    {
        private readonly IStiffnessDistributionOLD distribution;
        private readonly FetiDPDofSeparatorOLD dofSeparator;
        private readonly IModel model;

        public FetiDPSubdomainGlobalMappingOLD(IModel model, FetiDPDofSeparatorOLD dofSeparator,
            IStiffnessDistributionOLD distribution)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.distribution = distribution;
        }

        public double CalculateGlobalForcesNorm(Dictionary<int, IVectorView> subdomainForces)
        {
            //TODO: This can be optimized: calculate the dot product f*f for the internal dofs of each subdomain separately,
            //      only assemble global vector for the boundary dofs, find its dot product with itself, add the contributions
            //      for the internal dofs and finally apply SQRT(). This would greatly reduce the communication requirements.
            //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
            //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

            return GatherGlobalForces(subdomainForces).Norm2();
        }

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainRemainderDisplacements, 
            IVectorView globalCornerDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

            // Remainder dofs
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int id = subdomain.ID;
                int[] freeToGlobalDofs = model.GlobalDofOrdering.MapSubdomainToGlobalDofs(subdomain);
                int[] remainderToFreeDofs = dofSeparator.RemainderDofIndices[id];
                IVectorView remainderDisplacements = subdomainRemainderDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs are copied without averaging.
                foreach (int remainderDofIdx in dofSeparator.InternalDofIndices[id])
                {
                    int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
                    globalDisplacements[globalDofIdx] = remainderDisplacements[remainderDofIdx];
                }

                // For boundary dofs we take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int remainderDofIdx = dofSeparator.BoundaryDofIndices[id][i];
                    int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
                    globalDisplacements[globalDofIdx] += remainderDisplacements[remainderDofIdx] * boundaryDofCoeffs[i];
                }
            }

            // Corner dofs are copied without averaging.
            for (int i = 0; i < dofSeparator.NumGlobalCornerDofs; ++i)
            {
                int globalDofIdx = dofSeparator.GlobalCornerToFreeDofMap[i];
                globalDisplacements[globalDofIdx] = globalCornerDisplacements[i];
            }

            return globalDisplacements;
        }

        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
        {
            var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

            // Remainder dofs
            foreach (var subdomain in model.EnumerateSubdomains())
            {
                int id = subdomain.ID;
                int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapSubdomainToGlobalDofs(subdomain);
                int[] remainderToSubdomainDofs = dofSeparator.RemainderDofIndices[id];
                int[] cornerToSubdomainDofs = dofSeparator.CornerDofIndices[id];
                IVectorView freeDisplacements = subdomainDisplacements[id]; //TODO: benchmark the performance if this was concrete Vector

                // Internal dofs: We copy them without averaging.
                foreach (int remainderDofIdx in dofSeparator.InternalDofIndices[id])
                {
                    int subdomainDofIdx = remainderToSubdomainDofs[remainderDofIdx];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx];
                }

                // Boundary remainder dofs: We take average across subdomains with respect to multiplicity or stiffness. 
                double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
                for (int i = 0; i < dofSeparator.BoundaryDofIndices[id].Length; ++i)
                {
                    int subdomainDofIdx = remainderToSubdomainDofs[dofSeparator.BoundaryDofIndices[id][i]];
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                    globalDisplacements[globalDofIdx] += freeDisplacements[subdomainDofIdx] * boundaryDofCoeffs[i];
                }

                // Boundary corner dofs: We copy without averaging.
                foreach (int subdomainDofIdx in cornerToSubdomainDofs)
                {
                    int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];

                    // This will overwrite the value from a previous subdomain, but these values are the same.
                    globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx]; 
                }
            }

            return globalDisplacements;
        }

        public Vector GatherGlobalForces(Dictionary<int, IVectorView> subdomainForces)
        {
            var globalForces = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int id = subdomain.ID;
                int[] subdomainFreeToGlobalDofs = model.GlobalDofOrdering.MapSubdomainToGlobalDofs(subdomain);
                IVectorView forces = subdomainForces[id]; //TODO: benchmark the performance if this was concrete Vector

                for (int i = 0; i < forces.Length; ++i)
                {
                    // Internal forces will be copied (which is identical to adding 0 + single value).
                    // Boundary remainder forces will be summed. Previously we had distributed them depending on 
                    // homogeneity / heterogeneity (e.g. Ftot = 0.4 * Ftot + 0.6 * Ftot) and now we sum them. 
                    // Boundary corner forces are also summed. Previously we had also distributed them equally irregardless of 
                    // homogeneity / heterogeneity (e.g. Ftot = 0.5 * Ftot + 0.5 * Ftot) and now we sum them.
                    globalForces[subdomainFreeToGlobalDofs[i]] += forces[i];
                }
            }
            return globalForces;
        }
    }
}
