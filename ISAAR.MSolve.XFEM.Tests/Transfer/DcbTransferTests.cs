using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Entities;
using MPI;

namespace ISAAR.MSolve.XFEM.Tests.Transfer
{
    public static class DcbTransferTests
    {
        public const int numProcesses = 27;

        public static void TestEnrichmentTransfer(string[] args)
        {
            using (new MPI.Environment(ref args))
            {
                int master = 0;
                //int[] processesToSubdomains = Enumerable.Range(0, numProcesses).ToArray();
                int[][] processesToSubdomains = new int[numProcesses][];
                for (int p = 0; p < numProcesses; ++p) processesToSubdomains[p] = new int[] { p };
                var procs = new ProcessDistribution(Communicator.world, master, processesToSubdomains);

                // Expected model
                DcbBenchmarkBelytschko expectedBenchmark = CreateExpectedModel();
                XModel expectedModel = expectedBenchmark.Model;
                expectedModel.ConnectDataStructures();
                TrackingExteriorCrackLsm expectedCrack = expectedBenchmark.Crack;

                // Actual model
                DcbBenchmarkBelytschkoMpi actualBenchmark = CreateActualModel(procs);
                IXModelMpi actualModel = actualBenchmark.Model;
                actualModel.ConnectDataStructures();
                actualModel.ScatterSubdomains();
                var actualCrack = (TrackingExteriorCrackLsmMpiCentralized)actualBenchmark.Crack;
                //Console.WriteLine($"Process {procs.OwnRank}: Scattering lsm data");
                actualCrack.ScatterCrackData(actualModel);

                // Check that everything is in its correct subdomain and has the correct state.
                //Console.WriteLine($"Process {procs.OwnRank}: Checking lsm data for subdomain {procs.OwnSubdomainID}");
                int s = procs.GetSubdomainIDsOfProcess(procs.OwnRank).First();
                XSubdomain expectedSubdomain = expectedModel.Subdomains[s];
                SingleCrackLsm expectedLsm = expectedCrack.LevelSets;
                XSubdomain actualSubdomain = actualModel.GetXSubdomain(s);
                SingleCrackLsm actualLsm = actualCrack.LevelSets;
                EnrichmentComparions.CheckSameLevelSets(expectedLsm, expectedSubdomain, actualLsm, actualSubdomain);
                EnrichmentComparions.CheckSameEnrichments(expectedSubdomain, actualSubdomain);
            }
        }

        public static void TestModelTransfer(string[] args)
        {
            using (new MPI.Environment(ref args))
            {
                int master = 0;
                //int[] processesToSubdomains = Enumerable.Range(0, numProcesses).ToArray();
                int[][] processesToSubdomains = new int[numProcesses][];
                for (int p = 0; p < numProcesses; ++p) processesToSubdomains[p] = new int[] { p };
                var procs = new ProcessDistribution(Communicator.world, master, processesToSubdomains);

                // Expected model
                XModel expectedModel = CreateExpectedModel().Model;
                expectedModel.ConnectDataStructures();

                // Actual model
                IXModelMpi actualModel = CreateActualModel(procs).Model;
                actualModel.ConnectDataStructures();
                actualModel.ScatterSubdomains();

                // Check that all entities are in the correct subdomains and have the correct state.
                int s = procs.GetSubdomainIDsOfProcess(procs.OwnRank).First();
                XSubdomain actualSubdomain = actualModel.GetXSubdomain(s);
                XSubdomain expectedSubdomain = expectedModel.Subdomains[s];
                SubdomainComparisons.CheckSameNodes(expectedSubdomain, actualSubdomain);
                SubdomainComparisons.CheckSameElements(expectedSubdomain, actualSubdomain);
                SubdomainComparisons.CheckSameNodalLoads(expectedSubdomain, actualSubdomain);
                SubdomainComparisons.CheckSameNodalDisplacements(expectedSubdomain, actualSubdomain);
            }
        }

        private static DcbBenchmarkBelytschkoMpi CreateActualModel(ProcessDistribution procs)
        {
            int numElementsY = 15;
            int numSubdomainsY = 3;
            int numSubdomainsX = 3 * numSubdomainsY;

            var builder = new DcbBenchmarkBelytschkoMpi.Builder(procs, numElementsY, numSubdomainsX, numSubdomainsY);
            builder.HeavisideEnrichmentTolerance = 0.001;
            builder.MaxIterations = 8;
            builder.TipEnrichmentRadius = 0.0;
            builder.JintegralRadiusOverElementSize = 2.0;

            DcbBenchmarkBelytschkoMpi benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static DcbBenchmarkBelytschko CreateExpectedModel()
        {
            int numElementsY = 15;
            int numSubdomainsY = 3;
            int numSubdomainsX = 3 * numSubdomainsY;

            var builder = new DcbBenchmarkBelytschko.Builder(numElementsY, numSubdomainsX, numSubdomainsY);
            builder.HeavisideEnrichmentTolerance = 0.001;
            builder.MaxIterations = 8;
            builder.TipEnrichmentRadius = 0.0;
            builder.JintegralRadiusOverElementSize = 2.0;

            DcbBenchmarkBelytschko benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }
    }
}
