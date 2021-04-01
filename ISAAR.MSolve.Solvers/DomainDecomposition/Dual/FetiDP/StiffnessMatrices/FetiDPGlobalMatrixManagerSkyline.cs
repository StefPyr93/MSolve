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
    public class FetiDPGlobalMatrixManagerSkyline : FetiDPGlobalMatrixManagerBase
    {
        private readonly IReorderingAlgorithm reordering;
        private LdlSkyline inverseGlobalKccStar;

        public FetiDPGlobalMatrixManagerSkyline(IModel model, IFetiDPDofSeparator dofSeparator, IReorderingAlgorithm reordering):
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
            int[] skylineColHeights = FindSkylineColumnHeights(cornerNodeSelection);
            var skylineBuilder = SkylineBuilder.Create(dofSeparator.NumGlobalCornerDofs, skylineColHeights);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                IMatrixView subdomainKccStar = subdomainCoarseMatrices[subdomain];
                int[] subdomainToGlobalIndices = dofSeparator.GetCornerBooleanMatrix(subdomain).GetRowsToColumnsMap();
                skylineBuilder.AddSubmatrixSymmetric(subdomainKccStar, subdomainToGlobalIndices);
            }
            SkylineMatrix globalKccStar = skylineBuilder.BuildSkylineMatrix();

            // Factorization
            this.inverseGlobalKccStar = globalKccStar.FactorLdl(true);
        }

        protected override void ClearInverseCoarseProblemMatrixImpl() => inverseGlobalKccStar = null;

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector)
            => inverseGlobalKccStar.SolveLinearSystem(vector);

        private int[] FindSkylineColumnHeights(ICornerNodeSelection cornerNodeSelection)
        {
            // Only entries above the diagonal count towards the column height
            int[] colHeights = new int[dofSeparator.NumGlobalCornerDofs];
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> cornerNodes = cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain);

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = int.MaxValue;
                foreach (INode node in cornerNodes)
                {
                    foreach (int dof in dofSeparator.GlobalCornerDofOrdering.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (INode node in cornerNodes)
                {
                    foreach (int dof in dofSeparator.GlobalCornerDofOrdering.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return colHeights;
        }
    }
}
