using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.SchurComplements;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.SchurComplements
{
    public class SchurComplementCscTests
    {
        [Fact]
        public static void TestCalcComplementSymmetric()
        {
            // These are useful for debugging
            //string outputPath = @"C:\Users\Serafeim\Desktop\output.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var A = SkylineMatrix.CreateFromMatrix(Matrix.CreateFromArray(SymmPosDef10by10.Matrix));

            // Decompose the matrix into 4 parts. 
            //TODO: These should be hardcoded
            int[] indices2 = { 8, 3, 5, 1 }; // These will be condensed
            int[] indices1 = { 4, 0, 7, 2, 5, 6 };

            SkylineMatrix A22 = A.GetSubmatrixSymmetricSkyline(indices2);
            CscMatrix A21 = A.GetSubmatrixCsc(indices2, indices1);
            Matrix A11 = DenseStrategies.GetSubmatrix(A, indices1, indices1);
            SymmetricMatrix symA11 = A.GetSubmatrixSymmetricPacked(indices1);

            // Calculate the Schur complement
            LdlSkyline inverseA22 = A22.FactorLdl(false); //TODO: This should be hardcoded
            Matrix schurComplementExpected = SchurComplementCsc.CalcSchurComplementFull(A11, A21, inverseA22); //TODO: This should be hardcoded
            SymmetricMatrix schurComplementComputed1 = SchurComplementCsc.CalcSchurComplementSymmetric(symA11, A21, inverseA22);

            Matrix invA22TimesA21 = inverseA22.SolveLinearSystems(A21.CopyToFullMatrix());
            SymmetricMatrix schurComplementComputed2 = 
                SchurComplementCsc.CalcSchurComplementSymmetric(symA11, A21, invA22TimesA21);

            // Check
            Assert.True(schurComplementExpected.Equals(schurComplementComputed1));
            Assert.True(schurComplementExpected.Equals(schurComplementComputed2));
        }
    }
}
