using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed.Tests;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Performance;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using static ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests.PapagiannakisFetiDPTests2DMpi;

namespace ISAAR.MSolve.Solvers.Tests
{
    public static class Program
    {
        public static void Main(string[] args)
        {
            //ProfileFetiDPCantileverBeam2D.Run();

            var suite = new MpiTestSuite();
            RegisterFetiDP2dUnitTests(args, suite);
            //RegisterFetiDP2dIntegrationTests(args, suite);
            //RegisterFetiDP2dPapagiannakisTests(args, suite);
            suite.RunTests(args);
        }

        private static void RegisterFetiDP2dUnitTests(string[] args, MpiTestSuite suite)
        {
            //suite.AddFact(FetiDPDofSeparatorMpiTests.TestDofSeparation);
            //suite.AddFact(FetiDPDofSeparatorMpiTests.TestCornerBooleanMatrices);
            //suite.AddFact(FetiDPLagrangesEnumeratorMpiTests.TestBooleanMappingMatrices);
            //suite.AddFact(HomogeneousStiffnessDistributionMpiTests.TestBooleanMappingMatrices);

            //suite.AddFact(FetiDPMatrixManagerMpiTests.TestVectorsFbcFr);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKccKcrKrr, MatrixFormat.Dense);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKccKcrKrr, MatrixFormat.Skyline);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKccKcrKrr, MatrixFormat.SuiteSparse);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKbbKbiKii, MatrixFormat.Dense);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKbbKbiKii, MatrixFormat.Skyline);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestMatricesKbbKbiKii, MatrixFormat.SuiteSparse);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestStaticCondensations, MatrixFormat.Dense);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestStaticCondensations, MatrixFormat.Skyline);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestStaticCondensations, MatrixFormat.SuiteSparse);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestCoarseProblemMatrixAndRhs, MatrixFormat.Dense);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestCoarseProblemMatrixAndRhs, MatrixFormat.Skyline);
            //suite.AddTheory(FetiDPMatrixManagerMpiTests.TestCoarseProblemMatrixAndRhs, MatrixFormat.SuiteSparse);

            //suite.AddFact(FetiDPFlexibilityMatrixMpiTests.TestFIrcTimesVector);
            //suite.AddFact(FetiDPFlexibilityMatrixMpiTests.TestFIrcTransposedTimesVector);
            //suite.AddFact(FetiDPFlexibilityMatrixMpiTests.TestFIrrTimesVector);
            //suite.AddFact(FetiDPFlexibilityMatrixMpiTests.TestFIrrAndFIrcTransposedTimesVector);

            //suite.AddFact(FetiDPPreconditionerMpiTests.TestLumpedPreconditioner);
            //suite.AddFact(FetiDPPreconditionerMpiTests.TestDirichletPreconditioner);
            //suite.AddFact(FetiDPPreconditionerMpiTests.TestDiagonalDirichletPreconditioner);

            //suite.AddFact(FetiDPInterfaceProblemMpiTests.TestVectorDr);
            //suite.AddFact(FetiDPInterfaceProblemMpiTests.TestInterfaceProblemMatrix);
            //suite.AddFact(FetiDPInterfaceProblemMpiTests.TestInterfaceProblemRhs);
            //suite.AddFact(FetiDPInterfaceProblemMpiTests.TestInterfaceProblemSolution);

            suite.AddFact(FetiDPDisplacementsCalculatorMpiTests.TestCornerDisplacements);
            suite.AddFact(FetiDPDisplacementsCalculatorMpiTests.TestFreeDisplacements);

            //suite.AddFact(FetiDPSubdomainGlobalMappingMpiTests.TestGlobalDiplacements);
            //suite.AddFact(FetiDPSubdomainGlobalMappingMpiTests.TestGlobalForcesNorm);
        }

        private static void RegisterFetiDP2dIntegrationTests(string[] args, MpiTestSuite suite)
        {
            suite.AddTheory(FetiDPSolverMpiTests.TestSolutionSubdomainDisplacements, MatrixFormat.Skyline);
            suite.AddTheory(FetiDPSolverMpiTests.TestSolutionSubdomainDisplacements, MatrixFormat.SuiteSparse);
            suite.AddTheory(FetiDPSolverMpiTests.TestSolutionGlobalDisplacements, MatrixFormat.Skyline);
            suite.AddTheory(FetiDPSolverMpiTests.TestSolutionGlobalDisplacements, MatrixFormat.SuiteSparse);
        }

        private static void RegisterFetiDP2dPapagiannakisTests(string[] args, MpiTestSuite suite)
        {
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.Dirichlet, Residual.Approximate, 11, MatrixFormat.Skyline);
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.Dirichlet, Residual.Approximate, 11, MatrixFormat.SuiteSparse);
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.DirichletDiagonal, Residual.Approximate, 14, MatrixFormat.Skyline);
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.DirichletDiagonal, Residual.Approximate, 14, MatrixFormat.SuiteSparse);
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.Lumped, Residual.Approximate, 18, MatrixFormat.Skyline);
            suite.AddTheory(PapagiannakisFetiDPTests2DMpi.Run, 1.0, Precond.Lumped, Residual.Approximate, 18, MatrixFormat.SuiteSparse);
        }
    }
}
