using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: R1[s] = Br[s]^T * Q1[s]  is done in other classes as well. Perhaps I should do it once and store it instead of Q1.
//TODO: This only difference between this class and its 2D version is in the calculation of ur.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements
{
    public class FetiDP3dFreeDofDisplacementsCalculatorSerial : IFreeDofDisplacementsCalculator
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiDPMatrixManager matrixManager;
        private readonly IModel model;

        public FetiDP3dFreeDofDisplacementsCalculatorSerial(IModel model, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IAugmentationConstraints augmentationConstraints, 
            IFetiDPMatrixManager matrixManager)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.augmentationConstraints = augmentationConstraints;
            this.matrixManager = matrixManager;
        }

        public void CalculateSubdomainDisplacements(Vector lagranges, IFetiDPFlexibilityMatrix flexibility)
        {
            Vector ucTilde = CalcCornerDisplacementsAndAugmentedLagranges(flexibility, lagranges);
            Vector uc = ucTilde.GetSubvector(0, dofSeparator.NumGlobalCornerDofs);
            Vector mi = ucTilde.GetSubvector(dofSeparator.NumGlobalCornerDofs, ucTilde.Length);

            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                CalcAndStoreFreeDisplacements(subdomain, dofSeparator, lagrangesEnumerator, augmentationConstraints, 
                    matrixManager, lagranges, uc, mi);
            }
        }

        private Vector CalcCornerDisplacementsAndAugmentedLagranges(IFetiDPFlexibilityMatrix flexibility, Vector lagranges)
        {
            // ucTilde = inv(KccStarTilde) * (fcStarTilde + FIrcTilde^T * lagranges)
            Vector temp = flexibility.MultiplyFIrcTransposed(lagranges);
            temp.AddIntoThis(matrixManager.CoarseProblemRhs);
            return matrixManager.MultiplyInverseCoarseProblemMatrix(temp);
        }

        private static void CalcAndStoreFreeDisplacements(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IAugmentationConstraints augmentationConstraints,
            IFetiDPMatrixManager matrixManager, Vector lagranges, Vector uc, Vector mi)
        {
            IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Br = lagrangesEnumerator.GetBooleanMatrix(subdomain);
            IMappingMatrix Ba = augmentationConstraints.GetMatrixBa(subdomain);
            IMappingMatrix R1 = augmentationConstraints.GetMatrixR1(subdomain);

            // ur[s] = inv(Krr[s]) * (fr[s] - Br[s]^T * lambda - Krc[s] * Bc[s] * uc - Br[s]^T * Qr * mi) <=>
            // ur[s] = inv(Krr[s]) * (fr[s] - Br[s]^T * lambda - Krc[s] * Bc[s] * uc - R1[s] * Ba[s] * mi)
            Vector rhs = subdomainMatrices.Fr.Copy();
            rhs.SubtractIntoThis(Br.Multiply(lagranges, true));
            rhs.SubtractIntoThis(subdomainMatrices.MultiplyKrcTimes(Bc.Multiply(uc)));
            rhs.SubtractIntoThis(R1.Multiply(Ba.Multiply(mi)));
            Vector ur = subdomainMatrices.MultiplyInverseKrrTimes(rhs);

            // uf[s] = union(ur[s], ubc[s])
            // Remainder dofs
            var uf = Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs);
            int[] remainderDofs = dofSeparator.GetRemainderDofIndices(subdomain);
            uf.CopyNonContiguouslyFrom(remainderDofs, ur);

            // Corner dofs: ubc[s] = Bc[s] * uc
            Vector ubc = Bc.Multiply(uc);
            int[] cornerDofs = dofSeparator.GetCornerDofIndices(subdomain);
            uf.CopyNonContiguouslyFrom(cornerDofs, ubc);

            // Store uf[s]
            subdomainMatrices.LinearSystem.SolutionConcrete = uf;
        }
    }
}
