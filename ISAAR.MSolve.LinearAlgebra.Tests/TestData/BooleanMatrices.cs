using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// Sample matrices which only have +1, -1 or 0 entries.
    /// Authors: Serafeim Bakalakos 
    /// </summary>
    internal static class BooleanMatrices
    {
        internal static int[] MapCRowsToColumns => new int[] { 4, 0, 1, 8, 7 };


        internal static double[,] MatrixA10x5Single1PerCol => new double[,]
        {
            { 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0 },
            { 0, 0, 1, 0, 0 },
            { 0, 0, 0, 0, 1 },
            { 0, 0, 0, 0, 0 },
            { 1, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1, 0 }
        };

        internal static double[,] MatrixB5x10Single1PerColMultiplePerRow => new double[,]
        {
            { 0, 0, 0, 0, 1, 0, 1, 0, 0, 0 },
            { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
            { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1, 0, 1, 0, 0, 1, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }
        };

        internal static double[,] MatrixC5x10Single1PerRow => new double[,]
        {
            { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
            { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }
        };

        internal static double[,] MatrixD10x5Single1PerRowMultiplePerCol => new double[,]
        {
            { 1, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0 },
            { 1, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0 },
            { 0, 0, 0, 0, 1 },
            { 0, 1, 0, 0, 0 },
            { 0, 0, 1, 0, 0 },
            { 0, 0, 1, 0, 0 },
            { 0, 0, 0, 1, 0 },
            { 0, 0, 0, 1, 0 }
        };
    }
}
