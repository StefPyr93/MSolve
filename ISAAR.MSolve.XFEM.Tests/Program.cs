using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed.Tests;
using ISAAR.MSolve.XFEM.Tests.Transfer;

namespace ISAAR.MSolve.XFEM.Tests
{
    public static class Program
    {
        public static void Main(string[] args)
        {
            //Paper1.DoubleCantileverBeamRunner.RunTest();
            //Paper1.DoubleCantileverBeamRunner.RunVariousMeshes();
            //Paper1.DoubleCantileverBeamRunnerMpi.RunTest(args);

            //Paper1.HolesRunner.RunTest();
            Paper1.HolesRunnerMpi.RunTest(args);




            //DcbTransferTests.TestModelTransfer(args);
            //DcbTransferTests.TestEnrichmentTransfer(args);

            //HolesTransferTests.TestModelTransfer(args);
            //HolesTransferTests.TestEnrichmentTransfer(args);

            //RunMpiTests(args);

            //COMPDYN2019.DoubleCantileverBeam.Run();
            //COMPDYN2019.Holes.Run();
            //COMPDYN2019.Fillet.Run();
        }

        private static void RunMpiTests(string[] args)
        {
            var suite = new MpiTestSuite();

            //suite.AddFact(XModelMpiTests.TestModelTransfer, typeof(XModelMpiTests).Name, "TestModelTransfer");

            suite.RunTests(args);
        }
    }
}