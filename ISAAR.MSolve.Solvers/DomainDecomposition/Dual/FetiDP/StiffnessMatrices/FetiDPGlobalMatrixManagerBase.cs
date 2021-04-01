using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public abstract class FetiDPGlobalMatrixManagerBase : IFetiDPGlobalMatrixManager
    {
        protected readonly IFetiDPDofSeparator dofSeparator;
        protected readonly IModel model;

        private Vector globalFcStar;

        private bool hasInverseCoarseProblemMatrix;

        public FetiDPGlobalMatrixManagerBase(IModel model, IFetiDPDofSeparator dofSeparator)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
        }

        public Vector CoarseProblemRhs
        {
            get
            {
                if (globalFcStar == null) throw new InvalidOperationException( 
                    "The coarse problem RHS must be assembled from subdomains first.");
                return globalFcStar;
            }
        }

        public void CalcCoarseProblemRhs(Dictionary<ISubdomain, Vector> condensedRhsVectors) 
        {
            // globalFcStar = sum_over_s(Lc[s]^T * fcStar[s])
            globalFcStar = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain subdomain in condensedRhsVectors.Keys)
            {
                UnsignedBooleanMatrix Lc = dofSeparator.GetCornerBooleanMatrix(subdomain);
                Vector fcStar = condensedRhsVectors[subdomain];
                globalFcStar.AddIntoThisNonContiguouslyFrom(Lc.GetRowsToColumnsMap(), fcStar);
            }
        }

        public void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection, 
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices)
        {
            CalcInverseCoarseProblemMatrixImpl(cornerNodeSelection, subdomainCoarseMatrices);
            hasInverseCoarseProblemMatrix = true;
        }
        protected abstract void CalcInverseCoarseProblemMatrixImpl(ICornerNodeSelection cornerNodeSelection, 
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices);

        public void ClearCoarseProblemRhs() => globalFcStar = null;

        public void ClearInverseCoarseProblemMatrix()
        {
            ClearInverseCoarseProblemMatrixImpl();
            hasInverseCoarseProblemMatrix = false;
        }
        protected abstract void ClearInverseCoarseProblemMatrixImpl();

        public Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector)
        {
            if (!hasInverseCoarseProblemMatrix) throw new InvalidOperationException(
                "The inverse of the coarse problem matrix must be calculated first.");
            return MultiplyInverseCoarseProblemMatrixTimesImpl(vector);
        }
        protected abstract Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector);

        public abstract DofPermutation ReorderGlobalCornerDofs();
    }
}
