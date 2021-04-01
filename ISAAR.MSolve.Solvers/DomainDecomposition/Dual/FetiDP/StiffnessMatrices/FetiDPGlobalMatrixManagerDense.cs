using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPGlobalMatrixManagerDense : FetiDPGlobalMatrixManagerBase
    {
        private Matrix inverseGlobalKccStar;

        public FetiDPGlobalMatrixManagerDense(IModel model, IFetiDPDofSeparator dofSeparator) : base(model, dofSeparator)
        {
        }

        public override DofPermutation ReorderGlobalCornerDofs() => DofPermutation.CreateNoPermutation();

        protected override void CalcInverseCoarseProblemMatrixImpl(ICornerNodeSelection cornerNodeSelection,
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices)
        {
            // globalKccStar = sum_over_s(Bc[s]^T * KccStar[s] * Bc[s])
            var globalKccStar = Matrix.CreateZero(dofSeparator.NumGlobalCornerDofs, dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;
                IMatrixView subdomainKccStar = subdomainCoarseMatrices[subdomain];

                UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
                globalKccStar.AddIntoThis(Bc.ThisTransposeTimesOtherTimesThis(subdomainKccStar));
            }

            inverseGlobalKccStar = globalKccStar;
            inverseGlobalKccStar.InvertInPlace();
        }

        protected override void ClearInverseCoarseProblemMatrixImpl() => inverseGlobalKccStar = null;

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector) => inverseGlobalKccStar * vector;
    }
}
