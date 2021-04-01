using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This should not exist. Its code should be defined in FetiDPFlexibilityMatrixBase. It is not like other CPW Part classes
//      which actually stored state.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public class FetiDPSubdomainFlexibilityMatrix : IFetiDPSubdomainFlexibilityMatrix
    {
        //TODO:  If I store explicit matrices, then I would have to rebuild the flexibility matrix each time something changes. Not sure which is better
        //private readonly UnsignedBooleanMatrix Bc;
        //private readonly SignedBooleanMatrixColMajor Br;

        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly IFetiDPSubdomainMatrixManager matrixManager;
        private readonly ISubdomain subdomain;

        public FetiDPSubdomainFlexibilityMatrix(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator, IFetiDPMatrixManager matrixManager)
        {
            this.subdomain = subdomain;
            this.dofSeparator = dofSeparator;
            this.lagrangeEnumerator = lagrangeEnumerator;
            this.matrixManager = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
            //this.Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            //this.Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
        }

        public Vector MultiplyFIrc(Vector vector)
        {
            // FIrc[s] * x = sum_over_s( Br[s] * (inv(Krr[s]) * (Krc[s] * (Bc[s] * x))) ) 
            // Summing is delegated to another class.
            // This class performs: fIrc[s] * x = Br[s] * (inv(Krr[s]) * (Krc[s] * (Bc[s] * x)))

            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            #region debug1
            //Matrix Br_explicit = Br.MultiplyRight(Matrix.CreateIdentity(Br.NumColumns));
            // Bc_explicit = Bc.MultiplyRight(Matrix.CreateIdentity(Bc.NumColumns));
            var writer = new LinearAlgebra.Output.FullMatrixWriter();
            string pathBr = (new CnstValues()).debugString + @"\Inte_Test_Br_subd2.txt";
            string pathBc = (new CnstValues()).debugString + @"\Inte_Test_Bc_subd2.txt";
            //string pathBr = (new CnstValues()).debugString + @"\Main_Br_subd2.txt";
            //string pathBc = (new CnstValues()).debugString + @"\Main_Test_Bc_subd2.txt";

            //writer.WriteToFile(Br, pathBr);
            //writer.WriteToFile(Bc, pathBc);

            //Matrix Krc_explicit = MultiplyWithIdentity(Br.NumColumns, Bc.NumRows, (x, y) => y.CopyFrom(matrixManager.MultiplyKrcTimes(x)));
            // Bc_explicit = Bc.MultiplyRight(Matrix.CreateIdentity(Bc.NumColumns));
            #endregion
            Vector temp = Bc.Multiply(vector);
            temp = matrixManager.MultiplyKrcTimes(temp);
            temp = matrixManager.MultiplyInverseKrrTimes(temp);
            return Br.Multiply(temp);
        }

        public Vector MultiplyFIrcTransposed(Vector lagranges)
        {
            // FIrc[s]^T * x = sum_over_s( Bc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x))) ) 
            // Summing is delegated to another class.
            // This class performs: fIrc[s]^T * x = Bc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x)))

            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            Vector temp = Br.Multiply(lagranges, true);
            temp = matrixManager.MultiplyInverseKrrTimes(temp);
            temp = matrixManager.MultiplyKcrTimes(temp);
            return Bc.Multiply(temp, true);
        }

        public Vector MultiplyFIrr(Vector lagranges)
        {
            // FIrr[s] * x = sum_over_s( Br[s] * (inv(Krr[s]) * (Br[s]^T * x)) ) 
            // Summing is delegated to another class.
            // This class performs: fIrr[s] * x = Br[s] * (inv(Krr[s]) * (Br[s]^T * x))

            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            Vector temp = Br.Multiply(lagranges, true);
            temp = matrixManager.MultiplyInverseKrrTimes(temp);
            return Br.Multiply(temp);
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vector)
        {
            // Performs simultaneously a) fIrr[s] * x = Br[s] * (inv(Krr[s]) * (Br[s]^T * x))
            // and b) fIrc[s]^T * x = Bc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x)))
            // The computation of inv(Krr[s]) * (Br[s]^T * x) is common and can be reused.

            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Br = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            Vector invKrrTimesBrTimesVector = matrixManager.MultiplyInverseKrrTimes(Br.Multiply(vector, true));
            Vector FIrrTimesVector = Br.Multiply(invKrrTimesBrTimesVector);
            Vector FIrcTransposedTimesVector = Bc.Multiply(matrixManager.MultiplyKcrTimes(invKrrTimesBrTimesVector), true);
            return (FIrrTimesVector, FIrcTransposedTimesVector);
        }

        #region debug1 
        public static Matrix MultiplyWithIdentity(int numRows, int numCols, Action<Vector, Vector> matrixVectorMultiplication)
        {
            var result = Matrix.CreateZero(numRows, numCols);
            for (int j = 0; j < numCols; ++j)
            {
                var lhs = Vector.CreateZero(numCols);
                lhs[j] = 1.0;
                var rhs = Vector.CreateZero(numRows);
                matrixVectorMultiplication(lhs, rhs);
                result.SetSubcolumn(j, rhs);
            }
            return result;
        }
        #endregion
    }
}