using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs another version for boolean matrices, where the values array is omitted
//TODO: What happens if a column has only zero entries? Verify that the default 0 in rowIndices and values arrays good enough. 
//      Also add respective tests.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    /// <summary>
    /// Each column has exactly one non-zero entry.
    /// </summary>
    public class LocalToGlobalMappingMatrix : IMappingMatrix
    {
        /// <summary>
        /// The row index of the non-zero entry of each column
        /// </summary>
        private readonly int[] rowIndices;

        /// <summary>
        /// The value of the non-zero entry of each column
        /// </summary>
        private readonly double[] values;

        public LocalToGlobalMappingMatrix(int numRows, double[] values, int[] rowIndices)
        {
            this.rowIndices = rowIndices;
            this.values = values;
            this.NumRows = numRows;
            this.NumColumns = rowIndices.Length;
        }

        public int NumColumns { get; }

        public int NumRows { get; }

        public Matrix CopyToFullMatrix()
        {
            var result = new double[NumRows * NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[j * NumRows + rowIndices[j]] = values[j];
            }
            return Matrix.CreateFromArray(result, NumRows, NumColumns, false);
        }

        public void GetColumn(int j, Vector column)
        {
            column.Clear();
            column[rowIndices[j]] = values[j];
        }

        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyVectorTransposed(vector);
            else return MultiplyVectorUntransposed(vector);
        }

        public double MultiplyColumnTimes(Vector vector, int col)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            return values[col] * vector[rowIndices[col]];
        }

        public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyRightTransposed(other);
            else return MultiplyRightUntransposed(other);
        }

        /// <summary>
        /// Calculates the matrix product S = A21^T * inv(A22) * A21, where A21 = this, inv(A22) = other
        /// This method constructs inv(A22) * A21 one column at a time and uses that column to calculate the superdiagonal
        /// entries of the corresponding column of A21^T * inv(A22) * A21.
        /// </summary>
        public SymmetricMatrix MultiplyTransposeThisTimesOtherTimesThis(ITriangulation inverseA22)
        { //TODO: Unfortunately this cannot take advantage of MKL for CSC^T * vector.
            LocalToGlobalMappingMatrix A21 = this;
            var S = SymmetricMatrix.CreateZero(A21.NumColumns);

            var colA21 = Vector.CreateZero(A21.NumRows);
            for (int j = 0; j < A21.NumColumns; ++j)
            {
                // column j of (inv(A22) * A21) = inv(A22) * column j of A21
                A21.GetColumn(j, colA21);
                Vector colInvCA21 = inverseA22.SolveLinearSystem(colA21);

                // column j of (A21^T * inv(A22) * A21) = A21^T * column j of (inv(A22) * A21)
                // However we only need the superdiagonal part of this column. 
                // Thus we only multiply the rows i of A21^T (stored as columns i of A21) with i <= j. 
                for (int i = 0; i <= j; ++i)
                {
                    double dot = A21.MultiplyColumnTimes(colInvCA21, i);

                    // Assign the result for the current (i, j)
                    int indexS = S.Find1DIndex(i, j);
                    S.RawData[indexS] = dot;
                }
            }
            return S;
        }

        private Matrix MultiplyRightTransposed(Matrix other)
        {
            Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
            int numRowsResult = this.NumColumns;
            int numColsResult = other.NumColumns;
            var result = new double[numRowsResult * numColsResult];
            for (int col = 0; col < numColsResult; ++col)
            {
                int offset = col * numRowsResult;
                for (int j = 0; j < this.NumColumns; ++j)
                {
                    result[offset + j] = values[j] * other[rowIndices[j], col];
                }
            }
            return Matrix.CreateFromArray(result, numRowsResult, numColsResult, false);
        }

        private Matrix MultiplyRightUntransposed(Matrix other)
        {
            Preconditions.CheckMultiplicationDimensions(this.NumColumns, other.NumRows);
            int numRowsResult = this.NumRows;
            int numColsResult = other.NumColumns;
            var result = new double[numRowsResult * numColsResult];
            for (int col = 0; col < numColsResult; ++col)
            {
                int offset = col * numRowsResult;
                for (int j = 0; j < this.NumColumns; ++j)
                {
                    result[offset + rowIndices[j]] = values[j] * other[j, col];
                }
            }
            return Matrix.CreateFromArray(result, numRowsResult, numColsResult, false);
        }

        private Vector MultiplyVectorTransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[j] = values[j] * vector[rowIndices[j]];
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyVectorUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[rowIndices[j]] = values[j] * vector[j];
            }
            return Vector.CreateFromArray(result, false);
        }
    }
}
