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
    public class FetiDP3dGlobalMatrixManagerSuiteSparse : FetiDPGlobalMatrixManagerBase
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IReorderingAlgorithm reordering;

        private FetiDP3dCoarseProblemOrdering coarseOrdering;
        private CholeskySuiteSparse inverseGlobalKccStarTilde;

        public FetiDP3dGlobalMatrixManagerSuiteSparse(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints, IReorderingAlgorithm reordering)
            : base(model, dofSeparator)
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
            var dok = DokSymmetric.CreateEmpty(numCoarseDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int[] subToGlobalIndices = coarseOrdering.CoarseDofMapsSubdomainToGlobal[subdomain]; // subdomain-to-global mapping array
                IMatrixView subdomainKccStarTilde = subdomainCoarseMatrices[subdomain];
                dok.AddSubmatrixSymmetric(subdomainKccStarTilde, subToGlobalIndices);
            }
            SymmetricCscMatrix globalKccStar = dok.BuildSymmetricCscMatrix(true);

            // Factorization using SuiteSparse
            if (inverseGlobalKccStarTilde != null) inverseGlobalKccStarTilde.Dispose();
            inverseGlobalKccStarTilde = CholeskySuiteSparse.Factorize(globalKccStar, true);
        }

        protected override void ClearInverseCoarseProblemMatrixImpl()
        {
            if (inverseGlobalKccStarTilde != null) inverseGlobalKccStarTilde.Dispose();
            inverseGlobalKccStarTilde = null;
        }

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector)
        {
            //TODO: These should be handled by a dedicated PermutationMatrix class
            Vector permutedInputVector = vector.Reorder(coarseOrdering.CoarseDofPermutation, false);
            Vector permutedOutputVector = inverseGlobalKccStarTilde.SolveLinearSystem(permutedInputVector);
            Vector outputVector = permutedOutputVector.Reorder(coarseOrdering.CoarseDofPermutation, true);

            return outputVector;
        }
    }
}
