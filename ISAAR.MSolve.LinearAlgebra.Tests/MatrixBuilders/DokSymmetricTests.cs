using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.MatrixBuilders
{
    /// <summary>
    /// Tests for <see cref="DokSymmetric"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class DokSymmetricTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        private static int[] IndicesSet0 => new int[] { 0, 2, 4, 6, 12, 24, 32, 50, 64, 80 };
        private static int[] IndicesSet0Perm => new int[] { 32, 80, 64, 0, 12, 24, 6, 50, 4, 2 };
        private static int[] IndicesSet1 => new int[] { 4, 2, 5 };

        private static int[] IndicesSet2 => new int[] { 90, 10, 20, 60, 40, 50, 0, 70, 80, 30 };

        private static double[,] Matrix0 => MultiDiagonalMatrices.CreateSymmetric(100, new int[] { 2, 4, 8, 16, 32, 64 });

        private static double[,] Matrix1 => new double[,]
        {
            {  0,  0, 20,  0,  0,  0 },
            {  0,  1,  0, 31,  0,  0 },
            { 20,  0,  2,  0, 42,  0 },
            {  0, 31,  0,  3,  0, 53 },
            {  0,  0, 42,  0,  4,  0 },
            {  0,  0,  0, 53,  0,  5 }
        };
        

        private static DokSymmetric CreateDok(double[,] symmMatrix)
        {
            int n = symmMatrix.GetLength(0);
            var dok = DokSymmetric.CreateEmpty(n);
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (symmMatrix[i, j] != 0.0) dok[i, j] = symmMatrix[i, j];
                }
            }
            return dok;
        }

        [Fact]
        private static void TestAddSubmatrix()
        {
            var k1 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix1);
            var k2 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix2);
            var k3 = Matrix.CreateFromArray(GlobalMatrixAssembly.SubMatrix3);
            var expectedK = Matrix.CreateFromArray(GlobalMatrixAssembly.GlobalMatrix);

            var computedK = DokSymmetric.CreateEmpty(GlobalMatrixAssembly.GlobalOrder);
            computedK.AddSubmatrixSymmetric(k1, GlobalMatrixAssembly.IndicesDictionary1);
            computedK.AddSubmatrixSymmetric(k2, GlobalMatrixAssembly.IndicesDictionary2);
            computedK.AddSubmatrixSymmetric(k3, GlobalMatrixAssembly.IndicesDictionary3);

            comparer.AssertEqual(expectedK, computedK);
        }

        [Fact]
        private static void TestGetColumn()
        {
            Matrix dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            DokSymmetric dok = CreateDok(SparsePosDef10by10.Matrix);

            for (int j = 0; j < SparsePosDef10by10.Order; ++j)
            {
                comparer.AssertEqual(dense.GetColumn(j), dok.GetColumn(j)); //TODO: have hardcoded columns to compare against
            }
        }

        [Fact]
        private static void TestGetSubmatrixDokColMajor()
        {
            // These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;
            int[] rowIndices = indicesPerm;
            var colIndices = IndicesSet2;

            DokColMajor subMatrixCsc0 = matrixDok.GetSubmatrixDokColMajorNaive(indices, indices);
            DokColMajor subMatrixCsc1 = matrixDok.GetSubmatrixDokColMajor(indices, indices);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            
            Assert.True(subMatrixExpected.Equals(subMatrixCsc0));
            Assert.True(subMatrixExpected.Equals(subMatrixCsc1));

            DokColMajor subMatrixPermCsc0 = matrixDok.GetSubmatrixDokColMajorNaive(indicesPerm, indicesPerm);
            DokColMajor subMatrixPermCsc1 = matrixDok.GetSubmatrixDokColMajor(indicesPerm, indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermCsc0));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermCsc1));

            DokColMajor subMatrixRectCsc0 = matrixDok.GetSubmatrixDokColMajorNaive(rowIndices, colIndices);
            DokColMajor subMatrixRectCsc1 = matrixDok.GetSubmatrixDokColMajor(rowIndices, colIndices);
            Matrix subMatrixRectExpected = matrixFull.GetSubmatrix(rowIndices, colIndices);
            Assert.True(subMatrixRectExpected.Equals(subMatrixRectCsc0));
            Assert.True(subMatrixRectExpected.Equals(subMatrixRectCsc1));
        }

        [Fact]
        private static void TestGetSubmatrixFull()
        {
            // These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);
            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;
            int[] rowIndices = indicesPerm;
            var colIndices = IndicesSet2;

            Matrix subMatrixFull = matrixDok.GetSubmatrixFull(indices, indices);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            Assert.True(subMatrixExpected.Equals(subMatrixFull));

            Matrix subMatrixPermFull = matrixDok.GetSubmatrixFull(indicesPerm, indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermFull));

            Matrix subMatrixRectCscFull = matrixDok.GetSubmatrixFull(rowIndices, colIndices);
            Matrix subMatrixRectExpected = matrixFull.GetSubmatrix(rowIndices, colIndices);
            Assert.True(subMatrixRectExpected.Equals(subMatrixRectCscFull));
        }

        [Fact]
        private static void TestGetSubmatrixDok()
        {
            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;

            DokSymmetric subMatrixDok0 = matrixDok.GetSubmatrixSymmetricDokNaive(indices);
            DokSymmetric subMatrixDok1 = matrixDok.GetSubmatrixSymmetricDok(indices);
            //writer.WriteToFile(subMatrixSym, outputPath, true);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            Assert.True(subMatrixExpected.Equals(subMatrixDok0));
            Assert.True(subMatrixExpected.Equals(subMatrixDok1));

            DokSymmetric subMatrixPermDok0 = matrixDok.GetSubmatrixSymmetricDokNaive(indicesPerm);
            DokSymmetric subMatrixPermDok1 = matrixDok.GetSubmatrixSymmetricDok(indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermDok0));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermDok1));

            DokSymmetric matrix2 = CreateDok(new double[,]
            {
                {  0,  0, 20,  0,  0,  0 },
                {  0,  1,  0, 31,  0,  0 },
                { 20,  0,  2,  0, 42,  0 },
                {  0, 31,  0,  3,  0, 53 },
                {  0,  0, 42,  0,  4,  0 },
                {  0,  0,  0, 53,  0,  5 }
            });
            var rowsToKeep2 = new int[] { 4, 2, 5 };
            DokSymmetric submatrixExpected2 = CreateDok(new double[,]
            {
                {  4, 42, 0 }, 
                { 42,  2, 0 }, 
                {  0,  0, 5 }
            });

            DokSymmetric submatrixComputed2_0 = matrix2.GetSubmatrixSymmetricDokNaive(rowsToKeep2);
            DokSymmetric submatrixComputed2_1 = matrix2.GetSubmatrixSymmetricDok(rowsToKeep2);
            comparer.AssertEqual(submatrixExpected2, submatrixComputed2_0);
            comparer.AssertEqual(submatrixExpected2, submatrixComputed2_1);
        }

        [Fact]
        private static void TestGetSubmatrixSymmetricFull()
        {
            //// These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;

            Matrix subMatrixFull0 = matrixDok.GetSubmatrixSymmetricFullNaive(indices);
            Matrix subMatrixFull1 = matrixDok.GetSubmatrixSymmetricFull(indices);
            //writer.WriteToFile(subMatrixSym, outputPath, true);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            Assert.True(subMatrixExpected.Equals(subMatrixFull0));
            Assert.True(subMatrixExpected.Equals(subMatrixFull1));

            Matrix subMatrixPermFull0 = matrixDok.GetSubmatrixSymmetricFullNaive(indicesPerm);
            Matrix subMatrixPermFull1 = matrixDok.GetSubmatrixSymmetricFull(indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermFull0));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermFull1));
        }

        [Fact]
        private static void TestGetSubmatrixSymmetricPacked()
        {
            //// These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;

            SymmetricMatrix subMatrixPck0 = matrixDok.GetSubmatrixSymmetricPackedNaive(indices);
            SymmetricMatrix subMatrixPck1 = matrixDok.GetSubmatrixSymmetricPacked(indices);
            //writer.WriteToFile(subMatrixSym, outputPath, true);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            Assert.True(subMatrixExpected.Equals(subMatrixPck0));
            Assert.True(subMatrixExpected.Equals(subMatrixPck1));

            SymmetricMatrix subMatrixPermPck0 = matrixDok.GetSubmatrixSymmetricPackedNaive(indicesPerm);
            SymmetricMatrix subMatrixPermPck1 = matrixDok.GetSubmatrixSymmetricPacked(indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermPck0));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermPck1));
        }

        [Fact]
        private static void TestGetSubmatrixSymmetricPattern()
        {
            //// These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;

            SparsityPatternSymmetric subMatrixPattern0 = matrixDok.GetSubmatrixSymmetricPatternNaive(indices);
            SparsityPatternSymmetric subMatrixPattern1 = matrixDok.GetSubmatrixSymmetricPattern(indices);
            var subMatrixExpected = SparsityPatternSymmetric.CreateFromDense(matrixFull.GetSubmatrix(indices, indices));
            Assert.True(subMatrixExpected.Equals(subMatrixPattern0));
            Assert.True(subMatrixExpected.Equals(subMatrixPattern1));

            SparsityPatternSymmetric subMatrixPermPattern0 = matrixDok.GetSubmatrixSymmetricPatternNaive(indicesPerm);
            SparsityPatternSymmetric subMatrixPermPattern1 = matrixDok.GetSubmatrixSymmetricPattern(indicesPerm);
            var subMatrixPermExpected = 
                SparsityPatternSymmetric.CreateFromDense(matrixFull.GetSubmatrix(indicesPerm, indicesPerm));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermPattern0));
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermPattern1));
        }

        [Fact]
        private static void TestGetSubmatrixSymmetricSkyline()
        {
            //// These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var array2D = Matrix0;
            var matrixFull = Matrix.CreateFromArray(array2D);
            var matrixDok = DokSymmetric.CreateFromArray2D(array2D);

            int[] indices = IndicesSet0;
            int[] indicesPerm = IndicesSet0Perm;

            SkylineMatrix subMatrixSky = matrixDok.GetSubmatrixSymmetricSkyline(indices);
            //writer.WriteToFile(subMatrixSym, outputPath, true);
            Matrix subMatrixExpected = matrixFull.GetSubmatrix(indices, indices);
            Assert.True(subMatrixExpected.Equals(subMatrixSky));

            SkylineMatrix subMatrixPermSky = matrixDok.GetSubmatrixSymmetricSkyline(indicesPerm);
            Matrix subMatrixPermExpected = matrixFull.GetSubmatrix(indicesPerm, indicesPerm);
            Assert.True(subMatrixPermExpected.Equals(subMatrixPermSky));
        }

        [Fact] //TODO: Convert this to Theory
        private static void TestSplit_Full_DokColMajor()
        {
            var tries = new List<(double[,] matrix, int[] group0)>();
            tries.Add((Matrix0, IndicesSet0));
            tries.Add((Matrix0, IndicesSet0Perm));
            tries.Add((Matrix0, IndicesSet2));
            tries.Add((Matrix1, IndicesSet1));

            foreach ((double[,] matrixA, int[] group0) in tries)
            {
                var fullA = Matrix.CreateFromArray(matrixA);
                var dokA = DokSymmetric.CreateFromArray2D(matrixA);

                IEnumerable<int> allIndices = Enumerable.Range(0, fullA.NumColumns);
                int[] group1 = allIndices.Except(group0).ToArray();

                Matrix expectedA00 = fullA.GetSubmatrix(group0, group0);
                Matrix expectedA10 = fullA.GetSubmatrix(group1, group0);
                (Matrix A00, DokColMajor A10) = dokA.Split_Full_DokColMajor(group0, group1);

                comparer.AssertEqual(expectedA00, A00);
                comparer.AssertEqual(expectedA10, A10);


                Matrix expectedB00 = fullA.GetSubmatrix(group1, group1);
                Matrix expectedB10 = fullA.GetSubmatrix(group0, group1);
                (Matrix B00, DokColMajor B10) = dokA.Split_Full_DokColMajor(group1, group0);

                comparer.AssertEqual(expectedB00, B00);
                comparer.AssertEqual(expectedB10, B10);
            }
        }
        [Fact] //TODO: Convert this to Theory
        private static void TestSplit_Full_DokColMajor_DokSymmetric()
        {
            var tries = new List<(double[,] matrix, int[] group0)>();
            tries.Add((Matrix0, IndicesSet0));
            tries.Add((Matrix0, IndicesSet0Perm));
            tries.Add((Matrix0, IndicesSet2));
            tries.Add((Matrix1, IndicesSet1));

            foreach ((double[,] matrixA, int[] group0) in tries)
            {
                var fullA = Matrix.CreateFromArray(matrixA);
                var dokA = DokSymmetric.CreateFromArray2D(matrixA);

                IEnumerable<int> allIndices = Enumerable.Range(0, fullA.NumColumns);
                int[] group1 = allIndices.Except(group0).ToArray();

                Matrix expectedA00 = fullA.GetSubmatrix(group0, group0);
                Matrix expectedA10 = fullA.GetSubmatrix(group1, group0);
                Matrix expectedA11 = fullA.GetSubmatrix(group1, group1);
                (Matrix A00, DokColMajor A10, DokSymmetric A11) =
                    dokA.Split_Full_DokColMajor_DokSymmetric(group0, group1);

                comparer.AssertEqual(expectedA00, A00);
                comparer.AssertEqual(expectedA10, A10);
                comparer.AssertEqual(expectedA11, A11);


                Matrix expectedB00 = fullA.GetSubmatrix(group1, group1);
                Matrix expectedB10 = fullA.GetSubmatrix(group0, group1);
                Matrix expectedB11 = fullA.GetSubmatrix(group0, group0);
                (Matrix B00, DokColMajor B10, DokSymmetric B11) =
                    dokA.Split_Full_DokColMajor_DokSymmetric(group1, group0);

                comparer.AssertEqual(expectedB00, B00);
                comparer.AssertEqual(expectedB10, B10);
                comparer.AssertEqual(expectedB11, B11);
            }
        }

        [Fact] //TODO: Convert this to Theory
        private static void TestSplit_Packed_DokColMajor_DokSymmetric()
        {
            var tries = new List<(double[,] matrix, int[] group0)>();
            tries.Add((Matrix0, IndicesSet0));
            tries.Add((Matrix0, IndicesSet0Perm));
            tries.Add((Matrix0, IndicesSet2));
            tries.Add((Matrix1, IndicesSet1));

            foreach ((double[,] matrixA, int[] group0) in tries)
            {
                var fullA = Matrix.CreateFromArray(matrixA);
                var dokA = DokSymmetric.CreateFromArray2D(matrixA);

                IEnumerable<int> allIndices = Enumerable.Range(0, fullA.NumColumns);
                int[] group1 = allIndices.Except(group0).ToArray();

                Matrix expectedA00 = fullA.GetSubmatrix(group0, group0);
                Matrix expectedA10 = fullA.GetSubmatrix(group1, group0);
                Matrix expectedA11 = fullA.GetSubmatrix(group1, group1);
                (SymmetricMatrix A00, DokColMajor A10, DokSymmetric A11) =
                    dokA.Split_Packed_DokColMajor_DokSymmetric(group0, group1);
                
                comparer.AssertEqual(expectedA00, A00);
                comparer.AssertEqual(expectedA10, A10);
                comparer.AssertEqual(expectedA11, A11);


                Matrix expectedB00 = fullA.GetSubmatrix(group1, group1);
                Matrix expectedB10 = fullA.GetSubmatrix(group0, group1);
                Matrix expectedB11 = fullA.GetSubmatrix(group0, group0);
                (SymmetricMatrix B00, DokColMajor B10, DokSymmetric B11) =
                    dokA.Split_Packed_DokColMajor_DokSymmetric(group1, group0);

                comparer.AssertEqual(expectedB00, B00);
                comparer.AssertEqual(expectedB10, B10);
                comparer.AssertEqual(expectedB11, B11);
            }
        }

        [Fact]
        private static void TestIndexer()
        {
            Matrix dense = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            DokSymmetric dok = CreateDok(SparsePosDef10by10.Matrix);
            comparer.AssertEqual(dense, dok);
        }
    }
}
