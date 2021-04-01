using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="LocalToGlobalMappingMatrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LocalToGlobalMappingMatrixTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        internal static LocalToGlobalMappingMatrix CreateBooleanMatrix(double[,] matrix)
        {
            int m = matrix.GetLength(0);
            int n = matrix.GetLength(1);
            if (m < n) throw new ArgumentException("The matrix must have at least as many rows as columns");

            var values = new double[n];
            var rowIndices = new int[n];
            for (int j = 0; j < n; ++j)
            {
                bool isSet = false;
                for (int i = 0; i < m; ++i)
                {
                    if (matrix[i, j] != 0.0)
                    {
                        if (isSet) throw new ArgumentException("Only one entry per column can be non zero");
                        rowIndices[j] = i;
                        values[j] = matrix[i, j];
                        isSet = true;
                    }
                }
            }
            return new LocalToGlobalMappingMatrix(m, values, rowIndices);
        }

        [Fact]
        private static void TestMultiply()
        {
            var x = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(5));
            var y = Vector.CreateFromArray(RandomMatrices.CreateRandomVector(10));

            double[,] A = BooleanMatrices.MatrixA10x5Single1PerCol;
            LocalToGlobalMappingMatrix sparseΑ = CreateBooleanMatrix(A);
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
            var A = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(5, 12));
            var C = Matrix.CreateFromArray(RandomMatrices.CreateRandomMatrix(10, 4));

            double[,] B = BooleanMatrices.MatrixA10x5Single1PerCol;
            LocalToGlobalMappingMatrix sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            Matrix BA = sparseB.MultiplyRight(A);
            Matrix expectedBA = denseB.MultiplyRight(A);
            comparer.AssertEqual(expectedBA, BA);

            Matrix BtC = sparseB.MultiplyRight(C, true);
            Matrix expectedBtC = denseB.MultiplyRight(C, true);
            comparer.AssertEqual(expectedBtC, BtC);
        }

        [Fact]
        private static void TestMultiplyTransposeThisTimesOtherTimesThis()
        {
            var A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
            CholeskyFull factA = A.FactorCholesky(false);
            Matrix invA = A.Invert();

            double[,] B = BooleanMatrices.MatrixA10x5Single1PerCol;
            LocalToGlobalMappingMatrix sparseB = CreateBooleanMatrix(B);
            Matrix denseB = Matrix.CreateFromArray(B);

            SymmetricMatrix Bt_invA_B = sparseB.MultiplyTransposeThisTimesOtherTimesThis(factA);
            Matrix expectedBt_invA_B = denseB.MultiplyRight(invA, true) * denseB;
            comparer.AssertEqual(expectedBt_invA_B, Bt_invA_B);
        }
    }
}
