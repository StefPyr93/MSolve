using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Transfer;
using MPI;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests.Transfer
{
    public static class HolesTransferTests
    {
        public const int numProcesses = 10;

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
                HolesBenchmark expectedBenchmark = CreateExpectedModel();
                XModel expectedModel = expectedBenchmark.Model;
                expectedModel.ConnectDataStructures();

                // Actual model
                HolesBenchmarkMpi actualBenchmark = CreateActualModel(procs);
                IXModelMpi actualModel = actualBenchmark.Model;
                actualModel.ConnectDataStructures();
                actualModel.ScatterSubdomains();
                //Console.WriteLine($"Process {procs.OwnRank}: Scattering lsm data");
                actualBenchmark.Crack.ScatterCrackData(actualModel);

                // Check that everything is in its correct subdomain and has the correct state.
                //Console.WriteLine($"Process {procs.OwnRank}: Checking lsm data for subdomain {procs.OwnSubdomainID}");
                int s = procs.GetSubdomainIDsOfProcess(procs.OwnRank).First();
                XSubdomain expectedSubdomain = expectedModel.Subdomains[s];
                SingleCrackLsm expectedLeftLsm = expectedBenchmark.LeftCrack.LevelSets;
                SingleCrackLsm expectedRightLsm = expectedBenchmark.RightCrack.LevelSets;

                XSubdomain actualSubdomain = actualModel.GetXSubdomain(s);
                SingleCrackLsm actualLeftLsm = ((TrackingExteriorCrackLsmMpiCentralized)actualBenchmark.LeftCrack).LevelSets;
                SingleCrackLsm actualRightLsm = ((TrackingExteriorCrackLsmMpiCentralized)actualBenchmark.RightCrack).LevelSets;

                EnrichmentComparions.CheckSameLevelSets(expectedLeftLsm, expectedSubdomain, actualLeftLsm, actualSubdomain);
                EnrichmentComparions.CheckSameLevelSets(expectedRightLsm, expectedSubdomain, actualRightLsm, actualSubdomain);
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

        private static HolesBenchmarkMpi CreateActualModel(ProcessDistribution procs)
        {
            string meshPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.Parent.FullName
                + @"\Resources\holes_4442dofs.msh";

            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.
            var builder = new HolesBenchmarkMpi.Builder(procs, meshPath, growthLength);
            builder.HeavisideEnrichmentTolerance = 0.12;
            builder.MaxIterations = 12;
            builder.JintegralRadiusOverElementSize = 2.0;
            builder.TipEnrichmentRadius = 0.5;
            builder.BC = HolesBenchmarkMpi.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;
            HolesBenchmarkMpi benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel(numProcesses);
            return benchmark;
        }

        private static HolesBenchmark CreateExpectedModel()
        {
            string meshPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.Parent.FullName
                + @"\Resources\holes_4442dofs.msh";

            HolesBenchmark benchmark = HolesBenchmark.CreateMultiSubdomainBenchmark(10, () =>
            {
                double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.
                var builder = new HolesBenchmark.Builder(meshPath, growthLength);
                builder.HeavisideEnrichmentTolerance = 0.12;
                builder.MaxIterations = 12;
                builder.JintegralRadiusOverElementSize = 2.0;
                builder.TipEnrichmentRadius = 0.5;
                builder.BC = HolesBenchmark.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;
                HolesBenchmark singleSubdomainBenchmark = builder.BuildBenchmark();
                singleSubdomainBenchmark.InitializeModel();
                return singleSubdomainBenchmark;
            });

            return benchmark;
        }
    }
}
