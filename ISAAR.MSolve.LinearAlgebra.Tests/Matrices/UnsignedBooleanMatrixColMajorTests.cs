using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="UnsignedBooleanMatrixColMajor"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class UnsignedBooleanMatrixColMajorTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        internal static UnsignedBooleanMatrixColMajor CreateBooleanMatrix(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            var booleanMatrix = new UnsignedBooleanMatrixColMajor(m, n);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (matrix[i, j] == 1.0) booleanMatrix.AddEntry(i, j);
                    else if (matrix[i, j] != 0.0) throw new ArgumentException("Not a valid unsigned boolean matrix");
                }
            }
            return booleanMatrix;
        }

        [Fact]
        private static void TestGetColumns()
        {
            double[,] A = BooleanMatrices.MatrixD10x5Single1PerRowMultiplePerCol;
            UnsignedBooleanMatrixColMajor sparseA = CreateBooleanMatrix(A);
            var denseA = Matrix.CreateFromArray(A);

            int[] colsToKeepOption0 = { 1, 3 };
            int[] colsToKeepOption1 = { 2, 0, 1 };
            int[] allRows = Enumerable.Range(0, sparseA.NumRows).ToArray();

            UnsignedBooleanMatrixColMajor submatrixOption0 = sparseA.GetColumns(colsToKeepOption0, true);
            Matrix submatrixOption0Expected = denseA.GetSubmatrix(allRows, colsToKeepOption0);
            comparer.AssertEqual(submatrixOption0Expected, submatrixOption0);

            UnsignedBooleanMatrixColMajor submatrixOption1 = sparseA.GetColumns(colsToKeepOption1, false);
            Matrix submatrixOption1Expected = denseA.GetSubmatrix(allRows, colsToKeepOption1);
            comparer.AssertEqual(submatrixOption1Expected, submatrixOption1);
        }

        [Fact]
        private static void TestGetRowsToColumnsMap()
        {
            double[,] A = BooleanMatrices.MatrixC5x10Single1PerRow;
            UnsignedBooleanMatrixColMajor sparseA = CreateBooleanMatrix(A);

            int[] mapComputed = sparseA.GetRowsToColumnsMap();
            int[] mapExpected = BooleanMatrices.MapCRowsToColumns;

            comparer.AssertEqual(mapExpected, mapComputed);
        }

        [Fact]
        private static void TestMultiplicationThisTransposeTimesOtherTimesThis()
        {
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(5, 5));

            double[,] B = BooleanMatrices.MatrixC5x10Single1PerRow;
            UnsignedBooleanMatrixColMajor sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            Matrix BtAB = sparseB.ThisTransposeTimesOtherTimesThis(A);
            Matrix expectedBtAB = denseB.ThisTransposeTimesOtherTimesThis(A);

            comparer.AssertEqual(expectedBtAB, BtAB);
        }

        [Fact]
        private static void TestMultiply()
        {
            var x = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(5));
            var y = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(10));

            double[,] A = BooleanMatrices.MatrixA10x5Single1PerCol;
            UnsignedBooleanMatrixColMajor sparseΑ = CreateBooleanMatrix(A);
            Matrix denseΑ = Matrix.CreateFromArray(A);

            Vector Ax = sparseΑ.Multiply(x);
            Vector expectedAx= denseΑ.Multiply(x);
            comparer.AssertEqual(expectedAx, Ax);

            Vector Aty = sparseΑ.Multiply(y, true);
            Vector expectedAty = denseΑ.Multiply(y, true);
            comparer.AssertEqual(expectedAty, Aty);
        }

        [Fact]
        private static void TestMultiplyLeft()
        {
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(16, 10));

            double[,] B = BooleanMatrices.MatrixA10x5Single1PerCol;
            UnsignedBooleanMatrixColMajor sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            Matrix AB = sparseB.MultiplyLeft(A);
            Matrix expectedAB = denseB.MultiplyLeft(A);
            comparer.AssertEqual(expectedAB, AB);
        }

        [Fact]
        private static void TestMultiplyRight()
        {
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(5, 12));
            var C = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(10, 4));

            double[,] B = BooleanMatrices.MatrixA10x5Single1PerCol;
            UnsignedBooleanMatrixColMajor sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            Matrix BA = sparseB.MultiplyRight(A);
            Matrix expectedBA = denseB.MultiplyRight(A);
            comparer.AssertEqual(expectedBA, BA);

            Matrix BtC = sparseB.MultiplyRight(C, true);
            Matrix expectedBtC = denseB.MultiplyRight(C, true);
            comparer.AssertEqual(expectedBtC, BtC);
        }
    }
}
