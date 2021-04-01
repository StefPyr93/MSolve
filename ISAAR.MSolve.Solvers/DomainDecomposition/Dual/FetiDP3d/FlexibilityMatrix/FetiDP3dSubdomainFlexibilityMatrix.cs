using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;

//TODO: This should not exist. Its code should be defined in FetiDPFlexibilityMatrixBase. It is not like other CPW Part classes
//      which actually stored state.
//TODO: Instead of extracting subvectors, implement multiplications that operate on subvectors
//TODO: Perhaps Qr should not exist in all processes. In this case R1[s] should be used, which was already calculated during the 
//      coarse problem matrix assembly.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix
{
    public class FetiDP3dSubdomainFlexibilityMatrix : IFetiDPSubdomainFlexibilityMatrix
    {
        //TODO:  If I store explicit matrices, then I would have to rebuild the flexibility matrix each time something changes. Not sure which is better
        //private readonly UnsignedBooleanMatrix Bc;
        //private readonly SignedBooleanMatrixColMajor Br;

        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly IFetiDPSubdomainMatrixManager matrixManager;
        private readonly ISubdomain subdomain;

        public FetiDP3dSubdomainFlexibilityMatrix(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator, IAugmentationConstraints augmentationConstraints, 
            IFetiDPMatrixManager matrixManager)
        {
            this.subdomain = subdomain;
            this.dofSeparator = dofSeparator;
            this.lagrangeEnumerator = lagrangeEnumerator;
            this.augmentationConstraints = augmentationConstraints;
            this.matrixManager = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
            //this.Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            //this.Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
        }

        public Vector MultiplyFIrc(Vector vector)
        {
            // FIrc_tilde[s] * x = sum_over_s(fIrc_tilde[s] * x) 
            // Summing is delegated to another class.
            // This class performs: fIrc_tilde[s] * x = Br[s] * (inv(Krr[s]) * ((Krc[s] * (Bc[s] * xc) + (Br[s]^T * (Qr * xm) ))

            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            IMappingMatrix Qr = augmentationConstraints.MatrixGlobalQr;

            // Krc[s] * (Bc[s] * xc)
            Vector xc = vector.GetSubvector(0, dofSeparator.NumGlobalCornerDofs);
            Vector temp = Bc.Multiply(xc);
            temp = matrixManager.MultiplyKrcTimes(temp);

            // Br[s]^T * (Qr * xm)
            Vector xm = vector.GetSubvector(dofSeparator.NumGlobalCornerDofs, 
                dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints);
            Vector temp2 = Qr.Multiply(xm);
            temp2 = Br.Multiply(temp2, true);

            // (Krc[s] * (Bc[s] * xc) + (Br[s] ^ T * (Qr * xm)
            temp.AddIntoThis(temp2);

            // Br[s] * (inv(Krr[s]) * previousSum
            temp = matrixManager.MultiplyInverseKrrTimes(temp);
            return Br.Multiply(temp);
        }

        public Vector MultiplyFIrcTransposed(Vector vector)
        {
            // FIrc[s]_tilde^T * x = sum_over_s(fIrc[s]^T * x) 
            // Summing is delegated to another class.
            // This class performs: fIrc_tilde[s]^T * x = [ Bc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x))) ; 
            //      Qr^T * (Br[s] * (inv(Krr[s]) * (Br[s]^T * x))) ]

            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            IMappingMatrix Qr = augmentationConstraints.MatrixGlobalQr;

            Vector Br_invKrr_x = Br.Multiply(vector, true);
            Br_invKrr_x = matrixManager.MultiplyInverseKrrTimes(Br_invKrr_x);

            Vector temp1 = matrixManager.MultiplyKcrTimes(Br_invKrr_x);
            temp1 = Bc.Multiply(temp1, true);

            Vector temp2 = Br.Multiply(Br_invKrr_x);
            temp2 = Qr.Multiply(temp2, true);

            return temp1.Append(temp2);
        }
        
        public Vector MultiplyFIrr(Vector vector)
        {
            // FIrr[s] * x = sum_over_s( Br[s] * (inv(Krr[s]) * (Br[s]^T * x)) ) 
            // Summing is delegated to another class.
            // This class performs: fIrr[s] * x = Br[s] * (inv(Krr[s]) * (Br[s]^T * x))

            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            Vector temp = Br.Multiply(vector, true);
            temp = matrixManager.MultiplyInverseKrrTimes(temp);
            return Br.Multiply(temp);
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vector)
        {
            // Performs simultaneously a) fIrr[s] * x = Br[s] * (inv(Krr[s]) * (Br[s]^T * x))
            // and b) fIrc_tilde[s]^T * x = [ Bc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x))) ; 
            //      Qr^T * (Br[s] * (inv(Krr[s]) * (Br[s]^T * x))) ]
            // The computation of inv(Krr[s]) * (Br[s]^T * x) is common and can be reused.

            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            IMappingMatrix Qr = augmentationConstraints.MatrixGlobalQr;

            // Common part
            Vector invKrrTimesBrTimesVector = matrixManager.MultiplyInverseKrrTimes(Br.Multiply(vector, true));

            // FIrr part
            Vector FIrrTimesVector = Br.Multiply(invKrrTimesBrTimesVector);

            //FIrc part
            Vector temp1 = matrixManager.MultiplyKcrTimes(invKrrTimesBrTimesVector);
            temp1 = Bc.Multiply(temp1, true);
            Vector temp2 = Br.Multiply(invKrrTimesBrTimesVector);
            temp2 = Qr.Multiply(temp2, true);
            Vector FIrcTransposedTimesVector = temp1.Append(temp2);

            return (FIrrTimesVector, FIrcTransposedTimesVector);
        }
    }
}
