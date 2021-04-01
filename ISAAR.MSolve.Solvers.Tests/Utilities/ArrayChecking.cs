using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.Utilities
{
    internal static class ArrayChecking
    {
        internal static void CheckEqual(int[] expected, int[] computed)
        {
            Assert.Equal(expected.Length, computed.Length);
            for (int i = 0; i < expected.Length; ++i) Assert.Equal(expected[i], computed[i]);
        }

        internal static void CheckEqualMpi(ProcessDistribution procs, int[] expected, int[] computed)
        {
            if (expected.Length != computed.Length)
            {
                throw new ArgumentException($"Process {procs.OwnRank}: Computed array length = {computed.Length},"
                    + $" expected array length = {expected.Length}");
            }
            for (int i = 0; i < expected.Length; ++i)
            {
                if (expected[i] == computed[i]) continue;
                throw new ArgumentException($"Process {procs.OwnRank}: computed[{i}] = {computed[i]},"
                    + " expected[{i}] = {expected[i]}");
            }
        }
    }
}
