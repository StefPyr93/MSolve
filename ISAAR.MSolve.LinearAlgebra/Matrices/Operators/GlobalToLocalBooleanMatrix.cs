using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs another version for boolean matrices, where the values array is omitted
//TODO: What happens if a column has only zero entries? Verify that the default 0 in rowIndices and values arrays good enough. 
//      Also add respective tests.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    /// <summary>
    /// Each row has exactly one 1 entry.
    /// </summary>
    public class GlobalToLocalBooleanMatrix : IMappingMatrix
    {
        //TODO: Perhaps this should also contain the number of cols as its last entry, in order to make MPI communication faster
        /// <summary>
        /// The column index of the non-zero entry of each row
        /// </summary>
        private readonly int[] colIndices;

        public GlobalToLocalBooleanMatrix(int numCols, int[] colIndices)
        {
            this.colIndices = colIndices;
            this.NumRows = colIndices.Length;
            this.NumColumns = numCols;
        }

        public int NumColumns { get; }

        public int NumRows { get; }

        public Matrix CopyToFullMatrix()
        {
            var result = new double[NumRows * NumColumns];
            for (int i = 0; i < NumRows; ++i)
            {
                result[colIndices[i] * NumRows + i] = 1;
            }
            return Matrix.CreateFromArray(result, NumRows, NumColumns, false);
        }

        //TODO: Perhaps rename this as GetNonZerosColumnIndices()
        //TODO: This is risky as it exposes internal data.
        public int[] GetRowsToColumnsMap() => colIndices;

        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyVectorTransposed(vector);
            else return MultiplyVectorUntransposed(vector);
        }

        public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyRightTransposed(other);
            else return MultiplyRightUntransposed(other);
        }

        public Matrix MultiplyThisTransposeTimesOtherTimesThis(IMatrixView other) //TODO: this should be implemented for symmetric matrices
        {
            //TODO: Move this class to project Solvers where the following assumptions are always satisfied.
            //TODO: Otherwise, rename this method to specify for which instances it works correctly.

            // Rows of this matrix correspond to rows of matrix "other" (local). Columns of this matrix correspond to columns
            // of matrix "result" (global).
            Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
            var result = Matrix.CreateZero(this.NumColumns, this.NumColumns);
            int[] otherToResultMap = colIndices;

            for (int otherCol = 0; otherCol < other.NumColumns; ++otherCol)
            {
                int resultCol = otherToResultMap[otherCol];
                for (int otherRow = 0; otherRow < other.NumRows; ++otherRow)
                {
                    int resultRow = otherToResultMap[otherRow];
                    result[resultRow, resultCol] = other[otherRow, otherCol];
                }
            }
            return result;
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
                for (int i = 0; i < NumRows; ++i)
                {
                    result[offset + colIndices[i]] = other[i, col];
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
                for (int i = 0; i < NumRows; ++i)
                {
                    result[offset + i] = other[colIndices[i], col];
                }
            }
            return Matrix.CreateFromArray(result, numRowsResult, numColsResult, false);
        }

        private Vector MultiplyVectorTransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            for (int i = 0; i < NumRows; ++i)
            {
                result[colIndices[i]] = vector[i];
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyVectorUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            for (int i = 0; i < NumRows; ++i)
            {
                result[i] = vector[colIndices[i]];
            }
            return Vector.CreateFromArray(result, false);
        }
    }
}
