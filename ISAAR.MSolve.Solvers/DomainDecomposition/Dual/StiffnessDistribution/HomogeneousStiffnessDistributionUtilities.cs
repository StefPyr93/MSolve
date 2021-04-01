using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    internal static class HomogeneousStiffnessDistributionUtilities
    {
        internal static double[] CalcBoundaryDofInverseMultiplicities(ISubdomain subdomain, 
            (INode node, IDofType dofType)[] boundaryDofs)
        {
            int s = subdomain.ID;
            Debug.WriteLine($"{new StackTrace().GetFrame(1).GetMethod().ReflectedType.Name}:"
                + $" Calculating the inverse multiplicities of the boundary dofs of subdomain {s}");
            var inverseMultiplicities = new double[boundaryDofs.Length];
            for (int i = 0; i < boundaryDofs.Length; ++i)
            {
                inverseMultiplicities[i] = 1.0 / boundaryDofs[i].node.Multiplicity;
            }
            return inverseMultiplicities;
        }


        //TODO: This should be modified to CSR or CSC format and then benchmarked against the implicit alternative.
        /// <summary>
        /// Calculates the product Bpb = Bb * inv(Mb) explicitly, stores it and uses it for multiplications.
        /// </summary>
        internal class ScalingBooleanMatrixExplicit : IMappingMatrix
        {
            private readonly Matrix explicitBpb;

            internal ScalingBooleanMatrixExplicit(double[] inverseBoundaryDofMultiplicities,
                SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
            {
                var invMb = DiagonalMatrix.CreateFromArray(inverseBoundaryDofMultiplicities);
                this.explicitBpb = boundarySignedBooleanMatrix.MultiplyRight(invMb.CopyToFullMatrix());
            }

            public int NumColumns => explicitBpb.NumColumns;

            public int NumRows => explicitBpb.NumRows;
            public Matrix CopyToFullMatrix() => LinearAlgebra.Commons.DenseStrategies.CopyToFullMatrix(this);

            public Vector Multiply(Vector vector, bool transposeThis = false)
                => explicitBpb.Multiply(vector, transposeThis);

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
                => explicitBpb.MultiplyRight(other, transposeThis);
        }

        /// <summary>
        /// Stores the matrices Bb and inv(Mb). Matrix-vector and matrix-matrix multiplications with Bpb = Bb * inv(Mb) are
        /// performed implicitly, e.g. Bpb * x = Bb * (inv(Mb) * x).
        /// </summary>
        internal class ScalingBooleanMatrixImplicit : IMappingMatrix
        {
            /// <summary>
            /// Signed boolean matrix with only the boundary dofs of the subdomain as columns. 
            /// </summary>
            private readonly SignedBooleanMatrixColMajor Bb;

            /// <summary>
            /// Inverse of the diagonal matrix that stores the multiplicity of each boundary dof of the subdomain.
            /// </summary>
            private readonly DiagonalMatrix invMb;

            internal ScalingBooleanMatrixImplicit(double[] inverseBoundaryDofMultiplicities,
                SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
            {
                this.invMb = DiagonalMatrix.CreateFromArray(inverseBoundaryDofMultiplicities);
                this.Bb = boundarySignedBooleanMatrix;
            }

            public int NumColumns => invMb.NumColumns;

            public int NumRows => Bb.NumRows;

            public Matrix CopyToFullMatrix() => LinearAlgebra.Commons.DenseStrategies.CopyToFullMatrix(this);


            public Vector Multiply(Vector vector, bool transposeThis = false)
            {
                if (transposeThis)
                {
                    // Bpb^T * x = (Bb * inv(Mb))^T * x = inv(Mb)^T * Bb^T * x = inv(Mb) * (Bb^T * x);
                    return invMb * Bb.Multiply(vector, true);
                }
                else
                {
                    // Bpb * x = Bb * (inv(Mb) * x)
                    return Bb.Multiply(invMb * vector);
                }
            }

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
            {
                if (transposeThis)
                {
                    // Bpb^T * X = (Bb * inv(Mb))^T * X = inv(Mb)^T * Bb^T * X = inv(Mb) * (Bb^T * X);
                    return invMb * Bb.MultiplyRight(other, true);
                }
                else
                {
                    // Bpb * X = Bb * (inv(Mb) * X)
                    return Bb.MultiplyRight(invMb * other);
                }
            }
        }
    }
}