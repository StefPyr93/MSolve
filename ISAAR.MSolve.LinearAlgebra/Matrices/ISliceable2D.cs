﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Can return subvectors and submatices containing select entries.
    /// Authors: Serafeim Bakalakos
    /// </summary>/// <summary>
    public interface ISliceable2D: IIndexable2D
    {
        /// <summary>
        /// Returns a submatrix with the rows and columns of the original matrix that correspond to the provided index arrays. 
        /// The relative order of the rows and columns in the returned subvector can be different than the original one. It is 
        /// defined by the order of the index arrays <paramref name="rowIndices"/> and <paramref name="colIndices"/>
        /// respectively.
        /// </summary>
        /// <param name="rowIndices">The indices of the rows that will be returned in the submatrix. Constraints:
        ///     0 &lt;= <paramref name="rowIndices"/>[i] &lt; <see cref="IIndexable2D.NumRows"/>, for all i.</param>
        /// <param name="colIndices">The indices of the columns that will be returned in the submatrix.Constraints:
        ///     0 &lt;= <paramref name="colIndices"/>[j] &lt; <see cref="IIndexable2D.NumColumns"/>, for all j.</param>
        /// <returns></returns>
        Matrix Slice(int[] rowIndices, int[] colIndices);

        /// <summary>
        /// Returns a submatrix with entries (i, j) of the original matrix, such that 
        /// <paramref name="rowStartInclusive"/> &lt;= i &lt; <paramref name="rowEndExclusive"/> and
        /// <paramref name="colStartInclusive"/> &lt;= j &lt; <paramref name="colEndExclusive"/>.
        /// </summary>
        /// <param name="rowStartInclusive">The index of the first row that will be returned in the submatrix. Constraints: 
        ///     0 &lt;= <paramref name="rowStartInclusive"/> &lt;= <paramref name="rowEndExclusive"/>.</param>
        /// <param name="rowEndExclusive">The index immediately after the last row that will be returned in the 
        ///     submatrix. Constraints: <paramref name="rowStartInclusive"/> &lt;= <paramref name="rowEndExclusive"/> &lt; 
        ///     <see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="colStartInclusive">The index of the first column that will be returned in the submatrix. Constraints: 
        ///     0 &lt;= <paramref name="colStartInclusive"/> &lt;= <paramref name="colStartInclusive"/>.</param>
        /// <param name="colEndExclusive">The index immediately after the last column that will be returned in the submatrix.
        ///     Constraints: <paramref name="colStartInclusive"/> &lt;= <paramref name="colEndExclusive"/> &lt; 
        ///     <see cref="IIndexable2D.NumColumns"/>.</param>
        /// <returns></returns>
        Matrix Slice(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive);

        /// <summary>
        /// Returns the column of the original matrix with index = <paramref name="colIndex"/>.
        /// </summary>
        /// <param name="colIndex">The index of the column to return. Constraints: 
        ///     0 &lt;= <paramref name="colIndex"/> &lt; <see cref="IIndexable2D.NumColumns"/>.</param>
        /// <returns></returns>
        Vector SliceColumn(int colIndex);

        /// <summary>
        /// Returns the row of the original matrix with index = <paramref name="rowIndex"/>.
        /// </summary>
        /// <param name="rowIndex">The index of the row to return. Constraints: 
        ///     0 &lt;= <paramref name="rowIndex"/> &lt; <see cref="IIndexable2D.NumRows"/>.</param>
        /// <returns></returns>
        Vector SliceRow(int rowIndex);
    }
}
