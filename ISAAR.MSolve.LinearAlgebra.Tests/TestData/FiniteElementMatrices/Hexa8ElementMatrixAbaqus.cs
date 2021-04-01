﻿using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData.FiniteElementMatrices
{
    /// <summary>
    /// Stiffness matrix of hexahedral 8-node element taken from Abaqus. The element is a cube. Various boundary 
    /// conditions are possible.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class Hexa8ElementMatrixAbaqus
    {
        internal static double[,] UnconstrainedStiffness => new double[,]
        {
            { 81704.05983, -26362.17949, -26362.17949, 24866.45299, -560.8974359, -24118.58974, 9909.188034, 10657.05128, 10657.05128, 24866.45299, -24118.58974, -560.8974359, -27857.90598, 560.8974359, 560.8974359, -38327.99145, 26362.17949, -10657.05128, -36832.26496, 24118.58974, 24118.58974, -38327.99145, -10657.05128, 26362.17949 },
            { -26362.17949, 81704.05983, 26362.17949, 560.8974359, -27857.90598, -560.8974359, -10657.05128, -38327.99145, -26362.17949, -24118.58974, 24866.45299, 560.8974359, -560.8974359, 24866.45299, 24118.58974, 26362.17949, -38327.99145, 10657.05128, 24118.58974, -36832.26496, -24118.58974, 10657.05128, 9909.188034, -10657.05128 },
            { -26362.17949, 26362.17949, 81704.05983, -24118.58974, 560.8974359, 24866.45299, -10657.05128, -26362.17949, -38327.99145, 560.8974359, -560.8974359, -27857.90598, -560.8974359, 24118.58974, 24866.45299, 10657.05128, -10657.05128, 9909.188034, 24118.58974, -24118.58974, -36832.26496, 26362.17949, 10657.05128, -38327.99145 },
            { 24866.45299, 560.8974359, -24118.58974, 81704.05983, 26362.17949, -26362.17949, 24866.45299, 24118.58974, -560.8974359, 9909.188034, -10657.05128, 10657.05128, -38327.99145, -26362.17949, -10657.05128, -27857.90598, -560.8974359, 560.8974359, -38327.99145, 10657.05128, 26362.17949, -36832.26496, -24118.58974, 24118.58974 },
            { -560.8974359, -27857.90598, 560.8974359, 26362.17949, 81704.05983, -26362.17949, 24118.58974, 24866.45299, -560.8974359, 10657.05128, -38327.99145, 26362.17949, -26362.17949, -38327.99145, -10657.05128, 560.8974359, 24866.45299, -24118.58974, -10657.05128, 9909.188034, 10657.05128, -24118.58974, -36832.26496, 24118.58974 },
            { -24118.58974, -560.8974359, 24866.45299, -26362.17949, -26362.17949, 81704.05983, 560.8974359, 560.8974359, -27857.90598, -10657.05128, 26362.17949, -38327.99145, 10657.05128, 10657.05128, 9909.188034, -560.8974359, -24118.58974, 24866.45299, 26362.17949, -10657.05128, -38327.99145, 24118.58974, 24118.58974, -36832.26496 },
            { 9909.188034, -10657.05128, -10657.05128, 24866.45299, 24118.58974, 560.8974359, 81704.05983, 26362.17949, 26362.17949, 24866.45299, 560.8974359, 24118.58974, -36832.26496, -24118.58974, -24118.58974, -38327.99145, 10657.05128, -26362.17949, -27857.90598, -560.8974359, -560.8974359, -38327.99145, -26362.17949, 10657.05128 },
            { 10657.05128, -38327.99145, -26362.17949, 24118.58974, 24866.45299, 560.8974359, 26362.17949, 81704.05983, 26362.17949, -560.8974359, -27857.90598, -560.8974359, -24118.58974, -36832.26496, -24118.58974, -10657.05128, 9909.188034, -10657.05128, 560.8974359, 24866.45299, 24118.58974, -26362.17949, -38327.99145, 10657.05128 },
            { 10657.05128, -26362.17949, -38327.99145, -560.8974359, -560.8974359, -27857.90598, 26362.17949, 26362.17949, 81704.05983, 24118.58974, 560.8974359, 24866.45299, -24118.58974, -24118.58974, -36832.26496, -26362.17949, 10657.05128, -38327.99145, 560.8974359, 24118.58974, 24866.45299, -10657.05128, -10657.05128, 9909.188034 },
            { 24866.45299, -24118.58974, 560.8974359, 9909.188034, 10657.05128, -10657.05128, 24866.45299, -560.8974359, 24118.58974, 81704.05983, -26362.17949, 26362.17949, -38327.99145, -10657.05128, -26362.17949, -36832.26496, 24118.58974, -24118.58974, -38327.99145, 26362.17949, 10657.05128, -27857.90598, 560.8974359, -560.8974359 },
            { -24118.58974, 24866.45299, -560.8974359, -10657.05128, -38327.99145, 26362.17949, 560.8974359, -27857.90598, 560.8974359, -26362.17949, 81704.05983, -26362.17949, 10657.05128, 9909.188034, 10657.05128, 24118.58974, -36832.26496, 24118.58974, 26362.17949, -38327.99145, -10657.05128, -560.8974359, 24866.45299, -24118.58974 },
            { -560.8974359, 560.8974359, -27857.90598, 10657.05128, 26362.17949, -38327.99145, 24118.58974, -560.8974359, 24866.45299, 26362.17949, -26362.17949, 81704.05983, -26362.17949, -10657.05128, -38327.99145, -24118.58974, 24118.58974, -36832.26496, -10657.05128, 10657.05128, 9909.188034, 560.8974359, -24118.58974, 24866.45299 },
            { -27857.90598, -560.8974359, -560.8974359, -38327.99145, -26362.17949, 10657.05128, -36832.26496, -24118.58974, -24118.58974, -38327.99145, 10657.05128, -26362.17949, 81704.05983, 26362.17949, 26362.17949, 24866.45299, 560.8974359, 24118.58974, 9909.188034, -10657.05128, -10657.05128, 24866.45299, 24118.58974, 560.8974359 },
            { 560.8974359, 24866.45299, 24118.58974, -26362.17949, -38327.99145, 10657.05128, -24118.58974, -36832.26496, -24118.58974, -10657.05128, 9909.188034, -10657.05128, 26362.17949, 81704.05983, 26362.17949, -560.8974359, -27857.90598, -560.8974359, 10657.05128, -38327.99145, -26362.17949, 24118.58974, 24866.45299, 560.8974359 },
            { 560.8974359, 24118.58974, 24866.45299, -10657.05128, -10657.05128, 9909.188034, -24118.58974, -24118.58974, -36832.26496, -26362.17949, 10657.05128, -38327.99145, 26362.17949, 26362.17949, 81704.05983, 24118.58974, 560.8974359, 24866.45299, 10657.05128, -26362.17949, -38327.99145, -560.8974359, -560.8974359, -27857.90598 },
            { -38327.99145, 26362.17949, 10657.05128, -27857.90598, 560.8974359, -560.8974359, -38327.99145, -10657.05128, -26362.17949, -36832.26496, 24118.58974, -24118.58974, 24866.45299, -560.8974359, 24118.58974, 81704.05983, -26362.17949, 26362.17949, 24866.45299, -24118.58974, 560.8974359, 9909.188034, 10657.05128, -10657.05128 },
            { 26362.17949, -38327.99145, -10657.05128, -560.8974359, 24866.45299, -24118.58974, 10657.05128, 9909.188034, 10657.05128, 24118.58974, -36832.26496, 24118.58974, 560.8974359, -27857.90598, 560.8974359, -26362.17949, 81704.05983, -26362.17949, -24118.58974, 24866.45299, -560.8974359, -10657.05128, -38327.99145, 26362.17949 },
            { -10657.05128, 10657.05128, 9909.188034, 560.8974359, -24118.58974, 24866.45299, -26362.17949, -10657.05128, -38327.99145, -24118.58974, 24118.58974, -36832.26496, 24118.58974, -560.8974359, 24866.45299, 26362.17949, -26362.17949, 81704.05983, -560.8974359, 560.8974359, -27857.90598, 10657.05128, 26362.17949, -38327.99145 },
            { -36832.26496, 24118.58974, 24118.58974, -38327.99145, -10657.05128, 26362.17949, -27857.90598, 560.8974359, 560.8974359, -38327.99145, 26362.17949, -10657.05128, 9909.188034, 10657.05128, 10657.05128, 24866.45299, -24118.58974, -560.8974359, 81704.05983, -26362.17949, -26362.17949, 24866.45299, -560.8974359, -24118.58974 },
            { 24118.58974, -36832.26496, -24118.58974, 10657.05128, 9909.188034, -10657.05128, -560.8974359, 24866.45299, 24118.58974, 26362.17949, -38327.99145, 10657.05128, -10657.05128, -38327.99145, -26362.17949, -24118.58974, 24866.45299, 560.8974359, -26362.17949, 81704.05983, 26362.17949, 560.8974359, -27857.90598, -560.8974359 },
            { 24118.58974, -24118.58974, -36832.26496, 26362.17949, 10657.05128, -38327.99145, -560.8974359, 24118.58974, 24866.45299, 10657.05128, -10657.05128, 9909.188034, -10657.05128, -26362.17949, -38327.99145, 560.8974359, -560.8974359, -27857.90598, -26362.17949, 26362.17949, 81704.05983, -24118.58974, 560.8974359, 24866.45299 },
            { -38327.99145, 10657.05128, 26362.17949, -36832.26496, -24118.58974, 24118.58974, -38327.99145, -26362.17949, -10657.05128, -27857.90598, -560.8974359, 560.8974359, 24866.45299, 24118.58974, -560.8974359, 9909.188034, -10657.05128, 10657.05128, 24866.45299, 560.8974359, -24118.58974, 81704.05983, 26362.17949, -26362.17949 },
            { -10657.05128, 9909.188034, 10657.05128, -24118.58974, -36832.26496, 24118.58974, -26362.17949, -38327.99145, -10657.05128, 560.8974359, 24866.45299, -24118.58974, 24118.58974, 24866.45299, -560.8974359, 10657.05128, -38327.99145, 26362.17949, -560.8974359, -27857.90598, 560.8974359, 26362.17949, 81704.05983, -26362.17949 },
            { 26362.17949, -10657.05128, -38327.99145, 24118.58974, 24118.58974, -36832.26496, 10657.05128, 10657.05128, 9909.188034, -560.8974359, -24118.58974, 24866.45299, 560.8974359, 560.8974359, -27857.90598, -10657.05128, 26362.17949, -38327.99145, -24118.58974, -560.8974359, 24866.45299, -26362.17949, -26362.17949, 81704.05983 }
        };
    }
}
