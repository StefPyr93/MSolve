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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: remove duplication between this and the 2D version
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices
{
    public class FetiDP3dGlobalMatrixManagerDense : FetiDPGlobalMatrixManagerBase
    {
        private readonly IAugmentationConstraints augmentationConstraints;

        private Matrix inverseGlobalKccStarTilde;

        public FetiDP3dGlobalMatrixManagerDense(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints) : base(model, dofSeparator)
        {
            this.augmentationConstraints = augmentationConstraints;
        }

        public override DofPermutation ReorderGlobalCornerDofs() => DofPermutation.CreateNoPermutation();

        protected override void CalcInverseCoarseProblemMatrixImpl(ICornerNodeSelection cornerNodeSelection, 
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices)
        {
            // globalKccStar = sum_over_s(Bc[s]^T * KccStar[s] * Bc[s])
            // globalKacStar = sum_over_s(Ba[s]^T * KacStar[s] * Bc[s])
            // globalKaaStar = sum_over_s(Ba[s]^T * KaaStar[s] * Ba[s])

            int numCorners = dofSeparator.NumGlobalCornerDofs;
            int numAugmented = augmentationConstraints.NumGlobalAugmentationConstraints;
            var globalKccStar = Matrix.CreateZero(numCorners, numCorners);
            var globalKacStar = Matrix.CreateZero(numAugmented, numCorners);
            var globalKaaStar = Matrix.CreateZero(numAugmented, numAugmented);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;

                // Split the matrix containing all coarse dofs of this subdomain to corner and augmented dofs
                int countC = dofSeparator.GetCornerDofIndices(subdomain).Length;
                int countA = augmentationConstraints.GetNumAugmentationDofs(subdomain);
                Matrix subdomainMatrix = subdomainCoarseMatrices[subdomain].CopyToFullMatrix();
                Matrix KccStar = subdomainMatrix.GetSubmatrix(0, countC, 0, countC);
                Matrix KaaStar = subdomainMatrix.GetSubmatrix(countC, countC + countA, countC, countC + countA);
                Matrix KacStar = subdomainMatrix.GetSubmatrix(countC, countC + countA, 0, countC);

                UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
                GlobalToLocalBooleanMatrix Ba = augmentationConstraints.GetMatrixBa(subdomain);

                globalKccStar.AddIntoThis(Bc.ThisTransposeTimesOtherTimesThis(KccStar));
                globalKacStar.AddIntoThis(Ba.MultiplyRight(Bc.MultiplyLeft(KacStar), true));
                globalKaaStar.AddIntoThis(Ba.MultiplyThisTransposeTimesOtherTimesThis(KaaStar));
            }

            // KccTilde = [Kcc, Kac'; Kac Kaa];
            Matrix topRows = globalKccStar.AppendRight(globalKacStar.Transpose());
            Matrix bottomRows = globalKacStar.AppendRight(globalKaaStar);
            inverseGlobalKccStarTilde = topRows.AppendBottom(bottomRows);

            // Invert
            inverseGlobalKccStarTilde.InvertInPlace();
        }

        protected override void ClearInverseCoarseProblemMatrixImpl() => inverseGlobalKccStarTilde = null;

        protected override Vector MultiplyInverseCoarseProblemMatrixTimesImpl(Vector vector) 
            => inverseGlobalKccStarTilde * vector;
    }
}
