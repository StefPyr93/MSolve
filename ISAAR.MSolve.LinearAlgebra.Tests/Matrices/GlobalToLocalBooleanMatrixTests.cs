using System;
using System.Collections.Generic;
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
    /// Tests for <see cref="GlobalToLocalBooleanMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class GlobalToLocalBooleanMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        internal static GlobalToLocalBooleanMatrix CreateBooleanMatrix(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            if (m > n) throw new ArgumentException("The matrix must have at least as many columns as rows");

            var colIndices = new int[m];
            for (int i = 0; i < m; ++i)
            {
                bool isSet = false;
                for (int j = 0; j < n; ++j)
                {
                    if (matrix[i, j] != 0.0)
                    {
                        if (isSet) throw new ArgumentException("Only one entry per column can be non zero");
                        if (matrix[i, j] != 1) throw new ArgumentException("Only 0 or 1 entries are valid");
                        colIndices[i] = j;
                        isSet = true;
                    }
                }
            }
            return new GlobalToLocalBooleanMatrix(n, colIndices);
        }

        [Fact]
        private static void TestGetRowsToColumnsMap()
        {
            double[,] A = BooleanMatrices.MatrixC5x10Single1PerRow;
            GlobalToLocalBooleanMatrix sparseA = CreateBooleanMatrix(A);

            int[] mapComputed = sparseA.GetRowsToColumnsMap();
            int[] mapExpected = BooleanMatrices.MapCRowsToColumns;

            comparer.AssertEqual(mapExpected, mapComputed);
        }

        [Fact]
        private static void TestMultiplicationThisTransposeTimesOtherTimesThis()
        {
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(5, 5));

            double[,] B = BooleanMatrices.MatrixC5x10Single1PerRow;
            GlobalToLocalBooleanMatrix sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            Matrix BtAB = sparseB.MultiplyThisTransposeTimesOtherTimesThis(A);
            Matrix expectedBtAB = denseB.ThisTransposeTimesOtherTimesThis(A);

            comparer.AssertEqual(expectedBtAB, BtAB);
        }

        [Fact]
        private static void TestMultiply()
        {
            var x = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(10));
            var y = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(5));

            double[,] A = BooleanMatrices.MatrixC5x10Single1PerRow;
            GlobalToLocalBooleanMatrix sparseΑ = CreateBooleanMatrix(A);
            Matrix denseΑ = Matrix.CreateFromArray(A);

            Vector Ax = sparseΑ.Multiply(x);
            Vector expectedAx = denseΑ.Multiply(x);
            comparer.AssertEqual(expectedAx, Ax);

            Vector Aty = sparseΑ.Multiply(y, true);
            Vector expectedAty = denseΑ.Multiply(y, true);
            comparer.AssertEqual(expectedAty, Aty);
        }


        [Fact]
        private static void TestMultiplyRight()
        {
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(10, 4));
            var C = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(5, 12));

            double[,] B = BooleanMatrices.MatrixC5x10Single1PerRow;
            GlobalToLocalBooleanMatrix sparseB = CreateBooleanMatrix(B);
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
