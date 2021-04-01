using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;
using static ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.FetiDPMatrixManagerSerialTests;

//TODO: Remove duplication between this and the serial tests
//TODO: Mock all other FETI classes.
//TODO: Also check diagonal Kii. Actually there are a lot of missing stuff to check.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPMatrixManagerMpiTests
    {
        public static void TestCoarseProblemMatrixAndRhs(int numProcesses, MatrixFormat format)
        {
            IFetiDPMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) = 
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            FetiDPMatrixManagerMpi matrixManager = 
                PrepareCoarseProblemSubdomainMatrices(procs, model, dofSeparator, matricesFactory, format);

            // Calculate the global data to test
            matrixManager.CalcInverseCoarseProblemMatrix(Example4x4QuadsHomogeneous.DefineCornerNodeSelectionMpi(procs, model));
            matrixManager.CalcCoarseProblemRhs();

            if (procs.IsMasterProcess)
            {
                // Create explicit matrices from the matrix manager
                int numGlobalCornerDofs = dofSeparator.NumGlobalCornerDofs;
                Matrix globalInverseKccStar = ImplicitMatrixUtilities.MultiplyWithIdentity(numGlobalCornerDofs, numGlobalCornerDofs,
                    matrixManager.MultiplyInverseCoarseProblemMatrix);

                // Check
                double tol = 1E-13;
                Assert.True(Example4x4QuadsHomogeneous.MatrixGlobalKccStar.Invert().Equals(globalInverseKccStar, tol));
                Assert.True(Example4x4QuadsHomogeneous.VectorGlobalFcStar.Equals(matrixManager.CoarseProblemRhs, tol));
            }
        }

        public static void TestMatricesKbbKbiKii(int numProcesses, MatrixFormat format)
        {
            IFetiDPMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) =
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            var matrixManager = new FetiDPMatrixManagerMpi(procs, model, dofSeparator, matricesFactory);

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain sub = model.GetSubdomain(s);

                // Input data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                SetLinearSystemMatrix(subdomainMatrices.LinearSystem, Example4x4QuadsHomogeneous.GetMatrixKff(sub.ID), format);

                // Calculate the matrices to test
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractBoundaryInternalSubmatricesAndInvertKii(false);
                //subdomainMatrices.ExtractKbb();
                //subdomainMatrices.ExtractKbiKib();
                //subdomainMatrices.CalcInverseKii(false);

                // Create explicit matrices from the matrix manager
                int numBoundaryDofs = dofSeparator.GetBoundaryDofIndices(sub).Length;
                int numInternalDofs = dofSeparator.GetInternalDofIndices(sub).Length;
                Matrix Kbb = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numBoundaryDofs, numBoundaryDofs, subdomainMatrices.MultiplyKbbTimes);
                Matrix Kbi = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numBoundaryDofs, numInternalDofs, subdomainMatrices.MultiplyKbiTimes);
                Matrix inverseKii = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numInternalDofs, numInternalDofs, x => subdomainMatrices.MultiplyInverseKiiTimes(x, false));

                // Check
                double tol = 1E-13;
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKbb(sub.ID).Equals(Kbb, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKbi(sub.ID).Equals(Kbi, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKii(sub.ID).Invert().Equals(inverseKii, tol));
            }
        }

        public static void TestMatricesKccKcrKrr(int numProcesses, MatrixFormat format)
        {
            IFetiDPMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) =
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            var matrixManager = new FetiDPMatrixManagerMpi(procs, model, dofSeparator, matricesFactory);

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain sub = model.GetSubdomain(s);

                // Input data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                SetLinearSystemMatrix(subdomainMatrices.LinearSystem, Example4x4QuadsHomogeneous.GetMatrixKff(sub.ID), format);

                // Calculate the matrices to test
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.InvertKrr(true);

                // Create explicit matrices from the matrix manager
                int numCornerDofs = dofSeparator.GetCornerDofIndices(sub).Length;
                int numRemainderDofs = dofSeparator.GetRemainderDofIndices(sub).Length;
                Matrix Kcc = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numCornerDofs, numCornerDofs, subdomainMatrices.MultiplyKccTimes);
                Matrix Krc = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numRemainderDofs, numCornerDofs, subdomainMatrices.MultiplyKrcTimes);
                Matrix inverseKrr = ImplicitMatrixUtilities.MultiplyWithIdentity(
                    numRemainderDofs, numRemainderDofs, subdomainMatrices.MultiplyInverseKrrTimes);

                // Check
                double tol = 1E-13;
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKcc(sub.ID).Equals(Kcc, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKrc(sub.ID).Equals(Krc, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKrr(sub.ID).Invert().Equals(inverseKrr, tol));
            }
        }

        public static void TestStaticCondensations(int numProcesses, MatrixFormat format)
        {
            IFetiDPMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) =
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            var matrixManager = new FetiDPMatrixManagerMpi(procs, model, dofSeparator, matricesFactory);

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain sub = model.GetSubdomain(s);

                // Input data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                SetLinearSystemMatrix(subdomainMatrices.LinearSystem, Example4x4QuadsHomogeneous.GetMatrixKff(sub.ID), format);
                subdomainMatrices.LinearSystem.RhsConcrete = Example4x4QuadsHomogeneous.GetVectorFf(sub.ID);

                // Calculate the data to test
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
                subdomainMatrices.InvertKrr(true);
                subdomainMatrices.CalcCoarseProblemSubmatrices();
                subdomainMatrices.CalcCoarseProblemRhsSubvectors();

                // Check
                double tol = 1E-13;
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixKccStar(sub.ID).Equals(subdomainMatrices.CoarseProblemSubmatrix, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetVectorFcStar(sub.ID).Equals(subdomainMatrices.FcStar, tol));
            }
        }

        public static void TestVectorsFbcFr(int numProcesses)
        {
            IFetiDPMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(MatrixFormat.Dense);
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator) =
                FetiDPDofSeparatorMpiTests.CreateModelAndDofSeparator(numProcesses);

            var matrixManager = new FetiDPMatrixManagerMpi(procs, model, dofSeparator, matricesFactory);

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain sub = model.GetSubdomain(s);

                // Calculate the necessary vectors
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.LinearSystem.RhsConcrete = Example4x4QuadsHomogeneous.GetVectorFf(sub.ID);
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();

                // Check
                double tol = 1E-13;
                Assert.True(Example4x4QuadsHomogeneous.GetVectorFbc(sub.ID).Equals(subdomainMatrices.Fbc, tol));
                Assert.True(Example4x4QuadsHomogeneous.GetVectorFr(sub.ID).Equals(subdomainMatrices.Fr, tol));
            }
        }

        internal static FetiDPMatrixManagerMpi PrepareCoarseProblemSubdomainMatrices(ProcessDistribution procs, IModel model,
            IFetiDPDofSeparator dofSeparator, IFetiDPMatrixManagerFactory matricesFactory, MatrixFormat format)
        {
            var matrixManager = new FetiDPMatrixManagerMpi(procs, model, dofSeparator, matricesFactory);

            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain sub = model.GetSubdomain(s);
                // Input 
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                SetLinearSystemMatrix(subdomainMatrices.LinearSystem, Example4x4QuadsHomogeneous.GetMatrixKff(sub.ID), format);
                subdomainMatrices.LinearSystem.RhsConcrete = Example4x4QuadsHomogeneous.GetVectorFf(sub.ID);

                // Prepare the subdomain data
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
                subdomainMatrices.InvertKrr(true);
            }

            return matrixManager;
        }
    }
}
