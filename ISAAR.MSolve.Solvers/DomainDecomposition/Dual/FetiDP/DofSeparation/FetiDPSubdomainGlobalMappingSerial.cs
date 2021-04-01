using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPSubdomainGlobalMappingSerial
    {
        private readonly IStiffnessDistribution distribution;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly IModel model;

        public FetiDPSubdomainGlobalMappingSerial(IModel model, IFetiDPDofSeparator dofSeparator, 
            IStiffnessDistribution distribution)
        {
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

            return AssembleSubdomainVectors(getSubdomainForces).Norm2();
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

        //public Vector GatherGlobalDisplacementsOLD(Func<ISubdomain, Vector> getSubdomainRemainderDisplacements,
        //    IVectorView globalCornerDisplacements)
        //{
        //    var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

        //    // Remainder dofs
        //    foreach (ISubdomain subdomain in model.EnumerateSubdomains())
        //    {
        //        int[] freeToGlobalDofs = model.GlobalDofOrdering.MapSubdomainToGlobalDofs(subdomain);
        //        int[] remainderToFreeDofs = dofSeparator.GetRemainderDofIndices(subdomain);
        //        IVectorView remainderDisplacements = getSubdomainRemainderDisplacements(subdomain); //TODO: benchmark the performance if this was concrete Vector

        //        // Internal dofs are copied without averaging.
        //        foreach (int remainderDofIdx in dofSeparator.GetInternalDofIndices(subdomain))
        //        {
        //            int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
        //            globalDisplacements[globalDofIdx] = remainderDisplacements[remainderDofIdx];
        //        }

        //        // For boundary dofs we take average across subdomains with respect to multiplicity or stiffness. 
        //        double[] boundaryDofCoeffs = distribution.GetBoundaryDofCoefficients(subdomain);
        //        int[] boundaryDofIndices = dofSeparator.GetBoundaryDofIndices(subdomain);
        //        for (int i = 0; i < boundaryDofIndices.Length; ++i)
        //        {
        //            int remainderDofIdx = boundaryDofIndices[i];
        //            int globalDofIdx = freeToGlobalDofs[remainderToFreeDofs[remainderDofIdx]];
        //            globalDisplacements[globalDofIdx] += remainderDisplacements[remainderDofIdx] * boundaryDofCoeffs[i];
        //        }
        //    }

        //    // Corner dofs are copied without averaging.
        //    for (int i = 0; i < dofSeparator.NumGlobalCornerDofs; ++i)
        //    {
        //        int globalDofIdx = dofSeparator.GlobalCornerToFreeDofMap[i];
        //        globalDisplacements[globalDofIdx] = globalCornerDisplacements[i];
        //    }

        //    return globalDisplacements;
        //}

        //public Vector GatherGlobalDisplacementsOLD(Func<ISubdomain, Vector> getSubdomainFreeDisplacements)
        //{
        //    var globalDisplacements = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);

        //    // Remainder dofs
        //    foreach (var subdomain in model.EnumerateSubdomains())
        //    {
        //        int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapSubdomainToGlobalDofs(subdomain);
        //        int[] remainderToSubdomainDofs = dofSeparator.GetRemainderDofIndices(subdomain);
        //        int[] cornerToSubdomainDofs = dofSeparator.GetCornerDofIndices(subdomain);
        //        IVectorView freeDisplacements = getSubdomainFreeDisplacements(subdomain); //TODO: benchmark the performance if this was concrete Vector

        //        // Internal dofs: We copy them without averaging.
        //        foreach (int remainderDofIdx in dofSeparator.GetInternalDofIndices(subdomain))
        //        {
        //            int subdomainDofIdx = remainderToSubdomainDofs[remainderDofIdx];
        //            int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
        //            globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx];
        //        }

        //        // Boundary remainder dofs: We take average across subdomains with respect to multiplicity or stiffness. 
        //        double[] boundaryDofCoeffs = distribution.GetBoundaryDofCoefficients(subdomain);
        //        int[] boundaryDofIndices = dofSeparator.GetBoundaryDofIndices(subdomain);
        //        for (int i = 0; i < boundaryDofIndices.Length; ++i)
        //        {
        //            int subdomainDofIdx = remainderToSubdomainDofs[boundaryDofIndices[i]];
        //            int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
        //            globalDisplacements[globalDofIdx] += freeDisplacements[subdomainDofIdx] * boundaryDofCoeffs[i];
        //        }

        //        // Boundary corner dofs: We copy without averaging.
        //        foreach (int subdomainDofIdx in cornerToSubdomainDofs)
        //        {
        //            int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];

        //            // This will overwrite the value from a previous subdomain, but these values are the same.
        //            globalDisplacements[globalDofIdx] = freeDisplacements[subdomainDofIdx];
        //        }
        //    }

        //    return globalDisplacements;
        //}

        private Vector AssembleSubdomainVectors(Func<ISubdomain, Vector> getSubdomainForces)
        {
            var globalVector = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(subdomain, getSubdomainForces(subdomain), globalVector);
            }
            return globalVector;
        }
    }
}
