using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices
{
    public class FetiDP3dCoarseProblemOrdering
    {
        public FetiDP3dCoarseProblemOrdering(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints, IReorderingAlgorithm reordering)
        {
            // Find a naive global dof ordering of coarse dofs, where all corner dofs are numbered before all augmentation dofs
            CoarseDofOrdering = OrderCoarseProblemDofsNaively(dofSeparator, augmentationConstraints);

            // Reorder the coarse problem dofs
            (int[] permutation, bool oldToNew) = 
                ReorderCoarseProblemDofs(model, dofSeparator, augmentationConstraints, reordering, CoarseDofOrdering);
            CoarseDofOrdering.Reorder(permutation, oldToNew);
            CoarseDofPermutation = permutation;

            // Create subdomain to coarse problem (global) dof mappings so that they can be used (multiple times) 
            // during matrix assembly
            CoarseDofMapsSubdomainToGlobal = 
                MapCoarseDofsSubdomainToGlobal(model, dofSeparator, augmentationConstraints, CoarseDofOrdering);
        }

        public DofTable CoarseDofOrdering { get; }

        /// <summary>
        /// New-to-old
        /// </summary>
        public int[] CoarseDofPermutation { get; }

        public Dictionary<ISubdomain, int[]> CoarseDofMapsSubdomainToGlobal { get; }

        private static Dictionary<ISubdomain, int[]> MapCoarseDofsSubdomainToGlobal(IModel model, 
            IFetiDPDofSeparator dofSeparator, IAugmentationConstraints augmentationConstraints, DofTable coarseDofOrdering)
        {
            //TODO: This will be faster if I used the mappings defined by Bc, Ba and the coarse dof permutation array.
            var coarseDofMapsSubdomainToGlobal = new Dictionary<ISubdomain, int[]>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int numCornerDofs = dofSeparator.GetCornerDofIndices(subdomain).Length;
                int numAugmentationDofs = augmentationConstraints.GetNumAugmentationDofs(subdomain);
                var subToGlobalIndices = new int[numCornerDofs + numAugmentationDofs];
                coarseDofMapsSubdomainToGlobal[subdomain] = subToGlobalIndices;

                DofTable cornerDofOrdering = dofSeparator.GetCornerDofOrdering(subdomain);
                foreach ((INode node, IDofType dof, int localIdx) in cornerDofOrdering)
                {
                    subToGlobalIndices[localIdx] = coarseDofOrdering[node, dof];
                }

                DofTable augmentationDofOrdering = augmentationConstraints.GetAugmentationDofOrdering(subdomain);
                foreach ((INode node, IDofType dof, int localIdx) in augmentationDofOrdering)
                {
                    subToGlobalIndices[numCornerDofs + localIdx] = coarseDofOrdering[node, dof];
                }
            }
            return coarseDofMapsSubdomainToGlobal;
        }

        private static DofTable OrderCoarseProblemDofsNaively(IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints)
        {
            // Find a naive global dof ordering of coarse dofs, where all corner dofs are numbered before all augmentation dofs
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            DofTable coarseDofOrdering = dofSeparator.GlobalCornerDofOrdering.DeepCopy();
            foreach ((INode node, IDofType dof, int idx) in augmentationConstraints.GlobalAugmentationDofOrdering)
            { // TODO: Do this in DofTable
                coarseDofOrdering[node, dof] = numCornerDofs + idx;
            }
            return coarseDofOrdering;
        }

        private static (int[] permutation, bool oldToNew) ReorderCoarseProblemDofs(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints, IReorderingAlgorithm reordering, DofTable coarseDofOrdering)
        {
            // Find the sparsity pattern for this ordering
            int numCoarseProblemDofs =
                dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints;
            var pattern = SparsityPatternSymmetric.CreateEmpty(numCoarseProblemDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                // Treat each subdomain as a superelement with only its corner and midside nodes.
                int numSubdomainCornerDofs = dofSeparator.GetCornerDofIndices(subdomain).Length;
                int numSubdomainAugmentationDofs = augmentationConstraints.GetNumAugmentationDofs(subdomain);
                var subdomainToGlobalDofs = new int[numSubdomainCornerDofs + numSubdomainAugmentationDofs];

                // Corner dofs
                DofTable subdomainCornerDofOrdering = dofSeparator.GetCornerDofOrdering(subdomain);
                foreach ((INode node, IDofType dofType, int localIdx) in subdomainCornerDofOrdering)
                {
                    int globalIdx = coarseDofOrdering[node, dofType];
                    subdomainToGlobalDofs[localIdx] = globalIdx;
                }

                // Augmentation dofs follow corner dofs
                DofTable subdomainAugmentationDofOrdering = augmentationConstraints.GetAugmentationDofOrdering(subdomain);
                foreach ((INode node, IDofType dofType, int localIdx) in subdomainAugmentationDofOrdering)
                {
                    int globalIdx = coarseDofOrdering[node, dofType];
                    subdomainToGlobalDofs[localIdx + numSubdomainCornerDofs] = globalIdx;
                }

                pattern.ConnectIndices(subdomainToGlobalDofs, false);
            }

            // Reorder the coarse dofs
            (int[] permutation, bool oldToNew) result = reordering.FindPermutation(pattern);
            if (result.oldToNew) throw new NotImplementedException();
            return result;
        }
    }
}
