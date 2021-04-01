using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPGlobalMatrixManagerSuiteSparse : FetiDPGlobalMatrixManagerBase
    {
        private readonly IReorderingAlgorithm reordering;
        private CholeskySuiteSparse inverseGlobalKccStar;

        public FetiDPGlobalMatrixManagerSuiteSparse(IModel model, IFetiDPDofSeparator dofSeparator, 
            IReorderingAlgorithm reordering):
            base(model, dofSeparator)
        {
            this.reordering = reordering;
        }

        public override DofPermutation ReorderGlobalCornerDofs()
        {
            if (reordering == null) return DofPermutation.CreateNoPermutation();
            var pattern = SparsityPatternSymmetric.CreateEmpty(dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                // Treat each subdomain as a superelement with only its corner nodes.
                DofTable localCornerDofOrdering = dofSeparator.GetCornerDofOrdering(subdomain);
                int numLocalCornerDofs = localCornerDofOrdering.EntryCount;
                var subdomainToGlobalDofs = new int[numLocalCornerDofs];
                foreach ((INode node, IDofType dofType, int localIdx) in localCornerDofOrdering)
                {
                    int globalIdx = dofSeparator.GlobalCornerDofOrdering[node, dofType];
                    subdomainToGlobalDofs[localIdx] = globalIdx;
                }
                pattern.ConnectIndices(subdomainToGlobalDofs, false);
            }
            (int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
            return DofPermutation.Create(permutation, oldToNew);
        }

        protected override void CalcInverseCoarseProblemMatrixImpl(ICornerNodeSelection cornerNodeSelection,
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices)
        {
            // globalKccStar = sum_over_s(Lc[s]^T * KccStar[s] * Lc[s])

            // Assembly
            var dok = DokSymmetric.CreateEmpty(dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;
                IMatrixView subdomainKccStar = subdomainCoarseMatrices[subdomain];
                int[] subdomainToGlobalIndices = dofSeparator.GetCornerBooleanMatrix(subdomain).GetRowsToColumnsMap();
                dok.AddSubmatrixSymmetric(subdomainKccStar, subdomainToGlobalIndices);
            }
            SymmetricCscMatrix globalKccStar = dok.BuildSymmetricCscMatrix(true);

            // Factorization using SuiteSparse
            if (inverseGlobalKccStar != null) inverseGlobalKccStar.Dispose();
            inverseGlobalKccStar = CholeskySuiteSparse.Factorize(globalKccStar, true);
        }

        protected override void ClearInverseCoarseProblemMatrixImpl()
        {
            if (inverseGlobalKccStar != null) inverseGlobalKccStar.Dispose();
            inverseGlobalKccStar = null;
        }

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector)
            => inverseGlobalKccStar.SolveLinearSystem(vector);
    }
}
