using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: remove duplication between this and the 2D version. 
//TODO: remove duplication between Dense, Skyline and SuiteSparse versions.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices
{
    public class FetiDP3dGlobalMatrixManagerSkyline : FetiDPGlobalMatrixManagerBase
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IReorderingAlgorithm reordering;

        private FetiDP3dCoarseProblemOrdering coarseOrdering;
        private LdlSkyline inverseGlobalKccStarTilde;

        public FetiDP3dGlobalMatrixManagerSkyline(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints, IReorderingAlgorithm reordering) 
            : base (model, dofSeparator)
        {
            this.augmentationConstraints = augmentationConstraints;
            if (reordering == null) throw new NotImplementedException();
            this.reordering = reordering;
        }

        public override DofPermutation ReorderGlobalCornerDofs()
        {
            coarseOrdering = null;

            // For outside code the ordering remains the same
            return DofPermutation.CreateNoPermutation();
        }

        protected override void CalcInverseCoarseProblemMatrixImpl(ICornerNodeSelection cornerNodeSelection, 
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices)
        {
            // Order coarse dofs to minimize bandwidth
            //TODO: Needs thoughtful state management
            if (coarseOrdering == null)
            {
                coarseOrdering = new FetiDP3dCoarseProblemOrdering(model, dofSeparator, augmentationConstraints, reordering);
            }

            // globalKccStar = sum_over_s(Bc[s]^T * KccStar[s] * Bc[s])
            // globalKacStar = sum_over_s(Ba[s]^T * KacStar[s] * Bc[s])
            // globalKaaStar = sum_over_s(Ba[s]^T * KaaStar[s] * Ba[s])
            // However K matrices are joined as [KccStar, KacStar^T; KacStar, KaaStar] and Bc, Ba matrices are no longer very
            // useful, as the mappings they represent have been permuted.

            // Assembly
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numAugmentationDofs = augmentationConstraints.NumGlobalAugmentationConstraints;
            int numCoarseDofs = numCornerDofs + numAugmentationDofs;
            int[] skylineColHeights =
                FindSkylineColumnHeights(cornerNodeSelection, augmentationConstraints.MidsideNodesSelection);
            var skylineBuilder = SkylineBuilder.Create(numCoarseDofs, skylineColHeights);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int[] subToGlobalIndices = coarseOrdering.CoarseDofMapsSubdomainToGlobal[subdomain]; // subdomain-to-global mapping array
                IMatrixView subdomainMatrix = subdomainCoarseMatrices[subdomain]; // corner dofs followed by augmentation dofs
                skylineBuilder.AddSubmatrixSymmetric(subdomainMatrix, subToGlobalIndices);
            }
            SkylineMatrix globalKccStarTilde = skylineBuilder.BuildSkylineMatrix();

            // Factorization
            this.inverseGlobalKccStarTilde = globalKccStarTilde.FactorLdl(true);
        }

        protected override void ClearInverseCoarseProblemMatrixImpl() => inverseGlobalKccStarTilde = null;

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector)
        {
            Vector permutedInputVector = vector.Reorder(coarseOrdering.CoarseDofPermutation, false);
            Vector permutedOutputVector = inverseGlobalKccStarTilde.SolveLinearSystem(permutedInputVector);
            Vector outputVector = permutedOutputVector.Reorder(coarseOrdering.CoarseDofPermutation, true);

            return outputVector;
        }

        //TODO: Duplication between this and the 2D version
        private int[] FindSkylineColumnHeights(ICornerNodeSelection cornerNodeSelection, 
            IMidsideNodesSelection midsideNodesSelection)
        {
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            int numCoarseDofs = numCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints;

            // Only entries above the diagonal count towards the column height
            int[] colHeights = new int[numCoarseDofs];
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                // Put all coarse problem nodes (corners and midsides) together.
                var coarseNodes = new HashSet<INode>(cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain));
                coarseNodes.UnionWith(midsideNodesSelection.GetMidsideNodesOfSubdomain(subdomain));

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = int.MaxValue;
                foreach (INode node in coarseNodes)
                {
                    foreach (int dof in coarseOrdering.CoarseDofOrdering.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (INode node in coarseNodes)
                {
                    foreach (int dof in coarseOrdering.CoarseDofOrdering.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return colHeights;
        }
    }
}
