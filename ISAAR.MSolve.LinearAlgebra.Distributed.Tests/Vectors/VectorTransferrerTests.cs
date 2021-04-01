using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Tests;
using ISAAR.MSolve.LinearAlgebra.Distributed.Vectors;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;
using static ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer.TransferrerTestUtilities;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer
{
    public static class VectorTransferrerTests
    {
        private const int vectorLength = 10;

        /// <summary>
        /// All tests need 4 MPI processes.
        /// </summary>
        /// <param name="suite"></param>
        public static void RegisterAllTests(int numProcesses, MpiTestSuite suite)
        {
            // Tests for vector broadcasting
            suite.AddTheory(TestVectorBroadcast, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestVectorBroadcast, SubdomainDistribution.Uniform);
            suite.AddTheory(TestVectorBroadcast, SubdomainDistribution.Variable);

            // Tests for summing vectors
            suite.AddTheory(TestVectorSum, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestVectorSum, SubdomainDistribution.Uniform);
            suite.AddTheory(TestVectorSum, SubdomainDistribution.Variable);
            
            // Tests for summing vectors with lazy evaluation
            suite.AddTheory(TestVectorSumLazy, SubdomainDistribution.OnePerProcess);
            suite.AddTheory(TestVectorSumLazy, SubdomainDistribution.Uniform);
            suite.AddTheory(TestVectorSumLazy, SubdomainDistribution.Variable);
        }

        public static void TestVectorBroadcast(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            int s = 5;
            Vector vector = null;
            if (procs.IsMasterProcess) vector = GetSubdomainVector(s);
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref vector);
            
            if (!procs.IsMasterProcess)
            {
                Vector expected = GetSubdomainVector(s);
                Assert.True(expected.Equals(vector));
            }
        }

        public static void TestVectorSum(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            // Prepare vectors in each process
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            var processVectors = new Dictionary<int, Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                processVectors[s] = GetSubdomainVector(s);
            }

            // Sum the individual vectors
            var transferrer = new VectorTransferrer(procs);
            Vector sum1_master = transferrer.SumVectors(processVectors.Values);
            Vector sum2_master = null;
            if (procs.IsMasterProcess) sum2_master = Vector.CreateZero(vectorLength);
            transferrer.SumVectors(processVectors.Values, sum2_master);

            // Check
            if (procs.IsMasterProcess)
            {
                double tolerance = 1E-10;
                Vector sumExpected = GetTotalSum(procs);
                Assert.True(sumExpected.Equals(sum1_master, tolerance));
                Assert.True(sumExpected.Equals(sum2_master, tolerance));
            }
        }

        public static void TestVectorSumLazy(int numProcesses, SubdomainDistribution subdomainDistribution)
        {
            // Prepare vectors in each process to be lazily computed
            ProcessDistribution procs = DetermineProcesses(numProcesses, subdomainDistribution);
            IEnumerable<Vector> processVectors = procs.GetSubdomainIDsOfProcess(procs.OwnRank).Select(s => GetSubdomainVector(s));

            // Sum the individual vectors
            var transferrer = new VectorTransferrer(procs);
            Vector sum1_master = transferrer.SumVectors(processVectors);
            Vector sum2_master = null;
            if (procs.IsMasterProcess) sum2_master = Vector.CreateZero(vectorLength);
            transferrer.SumVectors(processVectors, sum2_master);

            // Check
            if (procs.IsMasterProcess)
            {
                double tolerance = 1E-10;
                Vector sumExpected = GetTotalSum(procs);
                Assert.True(sumExpected.Equals(sum1_master, tolerance));
                Assert.True(sumExpected.Equals(sum2_master, tolerance));
            }
        }

        private static Vector GetSubdomainVector(int subdomainID) => Vector.CreateWithValue(vectorLength, subdomainID);
        
        private static Vector GetTotalSum(ProcessDistribution procs)
        {
            double sum = 0.0;
            for (int p = 0; p < procs.Communicator.Size; ++p)
            {
                foreach (int s in procs.GetSubdomainIDsOfProcess(p)) sum += s;
            }
            return Vector.CreateWithValue(vectorLength, sum);
        }
    }
}
