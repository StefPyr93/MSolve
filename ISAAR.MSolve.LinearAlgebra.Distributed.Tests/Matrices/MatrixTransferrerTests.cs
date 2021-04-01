using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Tests;
using ISAAR.MSolve.LinearAlgebra.Distributed.Vectors;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;
using static ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer.TransferrerTestUtilities;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer
{
    public static class MatrixTransferrerTests
    {
        private const int numRows = 5, numCols = 4;

        /// <summary>
        /// All tests need 4 MPI processes.
        /// </summary>
        /// <param name="suite"></param>
        public static void RegisterAllTests(int numProcesses, MpiTestSuite suite)
        {
            // Tests for matrix broadcasting
            suite.AddTheory(TestMatrixBroadcast, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestMatrixBroadcast, SubdomainDistribution.Uniform);
            suite.AddTheory(TestMatrixBroadcast, SubdomainDistribution.Variable);

            // Tests for summing matrices
            suite.AddTheory(TestMatrixSum, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestMatrixSum, SubdomainDistribution.Uniform);
            suite.AddTheory(TestMatrixSum, SubdomainDistribution.Variable);

            // Tests for summing matrices with lazy evaluation
            suite.AddTheory(TestMatrixSumLazy, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestMatrixSumLazy, SubdomainDistribution.Uniform);
            suite.AddTheory(TestMatrixSumLazy, SubdomainDistribution.Variable);
        }

        public static void TestMatrixBroadcast(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            int s = 5;
            Matrix matrix = null;
            if (procs.IsMasterProcess) matrix = GetSubdomainMatrix(s);
            var transferrer = new MatrixTransferrer(procs);
            transferrer.BroadcastMatrix(ref matrix);
            
            if (!procs.IsMasterProcess)
            {
                Matrix expected = GetSubdomainMatrix(s);
                Assert.True(expected.Equals(matrix));
            }
        }

        public static void TestMatrixSum(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            // Prepare vectors in each process
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            var processMatrices = new Dictionary<int, Matrix>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                processMatrices[s] = GetSubdomainMatrix(s);
            }

            // Sum the individual vectors
            var transferrer = new MatrixTransferrer(procs);
            Matrix sum_master = null;
            if (procs.IsMasterProcess) sum_master = Matrix.CreateZero(numRows, numCols);
            transferrer.SumMatrices(processMatrices.Values, sum_master);

            // Check
            if (procs.IsMasterProcess)
            {
                double tolerance = 1E-10;
                Matrix sumExpected = GetTotalSum(procs);
                Assert.True(sumExpected.Equals(sum_master, tolerance));
                Assert.True(sumExpected.Equals(sum_master, tolerance));
            }
        }

        public static void TestMatrixSumLazy(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            // Prepare vectors in each process
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            IEnumerable<Matrix> processMatrices = procs.GetSubdomainIDsOfProcess(procs.OwnRank).Select(s => GetSubdomainMatrix(s));

            // Sum the individual vectors
            var transferrer = new MatrixTransferrer(procs);
            Matrix sum_master = null;
            if (procs.IsMasterProcess) sum_master = Matrix.CreateZero(numRows, numCols);
            transferrer.SumMatrices(processMatrices, sum_master);

            // Check
            if (procs.IsMasterProcess)
            {
                double tolerance = 1E-10;
                Matrix sumExpected = GetTotalSum(procs);
                Assert.True(sumExpected.Equals(sum_master, tolerance));
                Assert.True(sumExpected.Equals(sum_master, tolerance));
            }
        }

        private static Matrix GetSubdomainMatrix(int subdomainID) 
            => Matrix.CreateWithValue(numRows, numCols, subdomainID);
        
        private static Matrix GetTotalSum(ProcessDistribution procs)
        {
            double sum = 0.0;
            for (int p = 0; p < procs.Communicator.Size; ++p)
            {
                foreach (int s in procs.GetSubdomainIDsOfProcess(p)) sum += s;
            }
            return Matrix.CreateWithValue(numRows, numCols, sum);
        }
    }
}
