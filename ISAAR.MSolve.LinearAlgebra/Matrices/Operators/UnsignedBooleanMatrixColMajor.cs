using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    /// <summary>
    /// Sparse matrix with the non-zero entries being 1. Its main use is in domain decomposition solvers. 
    /// The internal data structures that store the non-zero entries are column major.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    [Serializable]
    public class UnsignedBooleanMatrixColMajor : IIndexable2D, IMappingMatrix
    {
        /// <summary>
        /// Non-zero entries: (row, columns). The key of the Dictionary is the column index, while the HashSet contains row 
        /// indices of non-zero entries.
        /// </summary>
        private readonly Dictionary<int, HashSet<int>> data;

        /// <summary>
        /// Initializes a new instance of <see cref="UnsignedBooleanMatrixColMajor"/> with the provided dimensions.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix. </param>
        public UnsignedBooleanMatrixColMajor(int numRows, int numColumns)
        {
            this.NumRows = numRows;
            this.NumColumns = numColumns;
            this.data = new Dictionary<int, HashSet<int>>(numColumns);
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        /// <remarks>
        /// The entries can be 0.0, 1.0 or -1.0
        /// </remarks>
        double IIndexable2D.this[int rowIdx, int colIdx] => this[rowIdx, colIdx];

        /// <summary>
        /// The entry with row index = rowIdx and column index = colIdx. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colIdx">The column index: 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">
        /// Thrown if <paramref name="rowIdx"/> or <paramref name="colIdx"/> violate the described constraints.
        /// </exception>
        public int this[int rowIdx, int colIdx]
        {
            get
            {
                if (data.TryGetValue(colIdx, out HashSet<int> rowIndices))
                {
                    if (rowIndices.Contains(rowIdx)) return 1;
                }
                return 0;
            }
        }

        /// <summary>
        /// Sets the entry with indices (<paramref name="colIdx"/>, <paramref name="rowIdx"/>) to 1.
        /// </summary>
        /// /// <param name="rowIdx">
        /// The column index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumColumns"/>.
        /// </param>
        /// <param name="colIdx">
        /// The row index of the entry to set. Constraints: 
        /// 0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumRows"/>.
        /// </param>
        public void AddEntry(int rowIdx, int colIdx)
        {
            if (data.TryGetValue(colIdx, out HashSet<int> rowIndices))
            {
                rowIndices.Add(rowIdx);
            }
            else
            {
                rowIndices = new HashSet<int>();
                data[colIdx] = rowIndices;
                rowIndices.Add(rowIdx);

            }
        }

        public Matrix CopyToFullMatrix() => DenseStrategies.CopyToFullMatrix((IMappingMatrix)this); //TODO: Optimize this

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public HashSet<int> GetNonZeroRowsOfColumn(int col) //TODO: This exposes internal data structures, which is risky
        {
            return data[col];
        }

        public UnsignedBooleanMatrixColMajor GetColumns(int[] colsToKeep, bool deepCopy) //TODO: This exposes internal data structures, which is risky
        {
            var clone = new UnsignedBooleanMatrixColMajor(this.NumRows, colsToKeep.Length);
            if (deepCopy)
            {
                for (int j = 0; j < colsToKeep.Length; ++j)
                {
                    var cloneColumn = new HashSet<int>();
                    foreach (int row in this.data[colsToKeep[j]]) cloneColumn.Add(row);
                    clone.data[j] = cloneColumn;
                }
            }
            else
            {
                for (int j = 0; j < colsToKeep.Length; ++j)
                {
                    clone.data[j] = this.data[colsToKeep[j]];
                }
            }
            return clone;
        }

        /// <summary>
        /// WARNING: this only works if this matrix is a mapping matrix L used in FETI solvers, meaning :
        /// 1) There are more columns than rows.
        /// 2) Each row of this matrix must have exactly one 1 and all other 0,
        /// 3) Each column must have at most one 1. It is possible that a column is completely 0.
        /// </summary>
        //TODO: This should be the actual way this matrix is stored. The dictionaries should be for building it only.
        //TODO: Perhaps rename this as GetNonZerosColumnIndices()
        public int[] GetRowsToColumnsMap()
        {
            int[] nonZeroRows = data.Keys.ToArray();
            foreach (var wholeCol in data)
            {
                Debug.Assert(wholeCol.Value.Count == 1);
                int col = wholeCol.Key;
                int row = wholeCol.Value.First();
                nonZeroRows[row] = col;
            }
            return nonZeroRows;
        }

        public Vector Multiply(Vector vector, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        public Matrix MultiplyLeft(IMatrixView other)
        {
            Preconditions.CheckMultiplicationDimensions(other, this);
            int numRowsResult = other.NumRows;
            int numColsResult = this.NumColumns;
            var result = new double[numRowsResult * numColsResult];
            foreach (var wholeCol in data)
            {
                int j = wholeCol.Key;
                int offset = j * numRowsResult;
                for (int i = 0; i < numRowsResult; ++i)
                {
                    double sum = 0.0;
                    foreach (int k in wholeCol.Value)
                    {
                        sum += other[i, k];
                    }
                    result[offset + i] = sum;
                }
            }
            return Matrix.CreateFromArray(result, numRowsResult, numColsResult, false);
        }

        public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
        {
            if (transposeThis) return MultiplyRightTransposed(other);
            else return MultiplyRightUntransposed(other);
        }

        /// <summary>
        /// WARNING: this only works if this matrix is a mapping matrix L used in FETI solvers, meaning :
        /// 1) There are more columns than rows.
        /// 2) Each row of this matrix must have exactly one 1 and all other 0,
        /// 3) Each column must have at most one 1. It is possible that a column is completely 0.
        /// Essentially this method performs global matrix assembly: globalK = L^T * localK * L.
        /// </summary>
        /// <param name="other"></param>
        public Matrix ThisTransposeTimesOtherTimesThis(IMatrixView other) //TODO: this should be implemented for symmetric matrices
        {
            //TODO: Move this class to project Solvers where the following assumptions are always satisfied.
            //TODO: Otherwise, rename this method to specify for which instances it works correctly.

            // Rows of this matrix correspond to rows of matrix "other" (local). Columns of this matrix correspond to columns
            // of matrix "result" (global).
            Preconditions.CheckMultiplicationDimensions(other.NumColumns, this.NumRows);
            var result = Matrix.CreateZero(this.NumColumns, this.NumColumns);
            int[] otherToResultMap = GetRowsToColumnsMap();

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

        private Vector MultiplyTransposed(Vector vector)
        {
            //TODO: I think that it will pay off to transpose an all integer CSR matrix and store both. Especially in the case 
            //     of subdomain boolean matrices, that little extra memory should not be of concern.
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            // Transpose it conceptually and multiply with the vector on the right. 
            foreach (var wholeCol in data)
            {
                double sum = 0.0;
                foreach (int row in wholeCol.Value)
                {
                    sum += vector[row];
                }
                result[wholeCol.Key] = sum;
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            foreach (var wholeCol in data)
            {
                double scalar = vector[wholeCol.Key];
                foreach (int row in wholeCol.Value)
                {
                    result[row] += scalar;
                }
            }
            return Vector.CreateFromArray(result, false);
        }

        private Matrix MultiplyRightTransposed(Matrix other)
        {
            //TODO: I think that it will pay off to transpose an all integer CSR matrix and store both. Especially in the case 
            //     of subdomain boolean matrices, that little extra memory should not be of concern.
            Preconditions.CheckMultiplicationDimensions(this.NumRows, other.NumRows);
            int numRowsResult = this.NumColumns;
            int numColsResult = other.NumColumns;
            var result = new double[numRowsResult * numColsResult];
            for (int j = 0; j < numColsResult; ++j)
            {
                int offset = j * numRowsResult;
                // Transpose it conceptually and multiply with the vector on the right. 
                foreach (var wholeCol in data)
                {
                    double sum = 0.0;
                    foreach (int row in wholeCol.Value)
                    {
                        sum += other[row, j];
                    }
                    result[offset + wholeCol.Key] = sum;
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
            for (int j = 0; j < numColsResult; ++j)
            {
                int offset = j * numRowsResult;
                foreach (var wholeCol in data)
                {
                    double scalar = other[wholeCol.Key, j];
                    foreach (int row in wholeCol.Value)
                    {
                        result[offset + row] += scalar;
                    }
                }
            }
            return Matrix.CreateFromArray(result, numRowsResult, other.NumColumns, false);
        }
    }
}
