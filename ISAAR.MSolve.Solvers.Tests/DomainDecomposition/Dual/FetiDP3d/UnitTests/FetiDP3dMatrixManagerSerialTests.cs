using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

//TODO: Mock all other FETI classes.
//TODO: Also check diagonal Kii. Actually there are a lot of missing stuff to check.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{ 
    public static class FetiDP3dMatrixManagerSerialTests
    {
        [Theory]
        [InlineData(MatrixFormat.Dense)]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestCoarseProblemMatrixAndRhs(MatrixFormat format)
        {
            // Create model
            Model model = ModelCreator.CreateModel();
            model.ConnectDataStructures();

            // Order free dofs
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Separate dofs and calculate the boolean matrices
            // Enumerate lagranges and calculate the boolean matrices

            ICornerNodeSelection cornerNodes = ModelCreator.DefineCornerNodeSelectionSerial(model);
            var dofSeparator = new FetiDPDofSeparatorSerial(model, cornerNodes);
            var crosspointStrategy = new FullyRedundantConstraints();
            var lagrangesEnumerator = new LagrangeMultipliersEnumeratorSerial(model, crosspointStrategy, dofSeparator);
            Dictionary<ISubdomain, HashSet<INode>> midsideNodes = ModelCreator.DefineMidsideNodesAll(model);
            IMidsideNodesSelection midsideNodesSelection = new UserDefinedMidsideNodes(midsideNodes,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
            IAugmentationConstraints augmentationConstraints =
                    new AugmentationConstraints(model, midsideNodesSelection, dofSeparator, lagrangesEnumerator);
            IFetiDP3dMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            var matrixManager = new FetiDP3dMatrixManagerSerial(model, dofSeparator, lagrangesEnumerator,
                augmentationConstraints, matricesFactory);


            SetKffRhs(model, matrixManager, format);
            //dofSeparator.SeparateDofs(new MockSeparatedDofReordering()); //TODO: This does not work as intended
            dofSeparator.SeparateDofs(matrixManager);
            lagrangesEnumerator.CalcBooleanMatrices(dofSeparator.GetRemainderDofOrdering);
            augmentationConstraints.CalcAugmentationMappingMatrices();
            PrepareCoarseProblemSubdomainMatrices(model, matrixManager, format);

            // Calculate the global data to test
            matrixManager.CalcInverseCoarseProblemMatrix(ModelCreator.DefineCornerNodeSelectionSerial(model));
            matrixManager.CalcCoarseProblemRhs();

            // Create explicit matrices from the matrix manager
            int coarseProblemSize = dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints;
            Matrix globalInverseKccStarTilde = ImplicitMatrixUtilities.MultiplyWithIdentity(coarseProblemSize,
                coarseProblemSize, matrixManager.MultiplyInverseCoarseProblemMatrix);

            // Check
            double tol = 1E-3;
            Assert.True(ExpectedGlobalMatrices.MatrixGlobalKccStarTildeSimple.Invert().Equals(globalInverseKccStarTilde, tol));
            Assert.True(ExpectedGlobalMatrices.VectorGlobalFcStarTildeSimple.Equals(matrixManager.CoarseProblemRhs, tol));
        }


        [Theory]
        [InlineData(MatrixFormat.Dense)]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestMatricesKbbKbiKii(MatrixFormat format)
        {
            IFetiDP3dMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            FetiDP3dMatrixManagerSerial matrixManager = PrepareCoarseProblemSubdomainMatrices(model, dofSeparator,
                lagrangesEnumerator, augmentationConstraints, matricesFactory, format);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Calculate the matrices to test
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractBoundaryInternalSubmatricesAndInvertKii(false);

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
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKbb(sub.ID).Equals(Kbb, tol));
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKbi(sub.ID).Equals(Kbi, tol));
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKii(sub.ID).Invert().Equals(inverseKii, tol));
            }
        }

        [Theory]
        [InlineData(MatrixFormat.Dense)]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestMatricesKcrKrr(MatrixFormat format)
        {
            IFetiDP3dMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            FetiDP3dMatrixManagerSerial matrixManager = PrepareCoarseProblemSubdomainMatrices(model, dofSeparator,
                lagrangesEnumerator, augmentationConstraints, matricesFactory, format);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Calculate the matrices to test
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
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
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKcc(sub.ID).Equals(Kcc, tol));
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKrc(sub.ID).Equals(Krc, tol));
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKrr(sub.ID).Invert().Equals(inverseKrr, tol));
            }
        }

        [Theory]
        [InlineData(MatrixFormat.Dense)]
        [InlineData(MatrixFormat.Skyline)]
        [InlineData(MatrixFormat.SuiteSparse)]
        public static void TestStaticCondensations(MatrixFormat format)
        {
            IFetiDP3dMatrixManagerFactory matricesFactory = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            FetiDP3dMatrixManagerSerial matrixManager = PrepareCoarseProblemSubdomainMatrices(model, dofSeparator,
                lagrangesEnumerator, augmentationConstraints, matricesFactory, format);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Calculate the data to test
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
                subdomainMatrices.InvertKrr(true);
                subdomainMatrices.CalcCoarseProblemSubmatrices();
                subdomainMatrices.CalcCoarseProblemRhsSubvectors();

                // Isolate KccStar
                Matrix KccStarTilde = subdomainMatrices.CoarseProblemSubmatrix.CopyToFullMatrix();
                int numCornerDofs = dofSeparator.GetCornerDofIndices(sub).Length;
                Matrix KccStar = KccStarTilde.GetSubmatrix(0, numCornerDofs, 0, numCornerDofs);

                // Check
                double tol = 1E-5;
                Assert.True(ExpectedSubdomainMatrices.GetMatrixKccStar(sub.ID).Equals(KccStar, tol));
                Assert.True(ExpectedSubdomainMatrices.GetVectorFcStar(sub.ID).Equals(subdomainMatrices.FcStar, tol));
            }
        }

        [Fact]
        public static void TestVectorDr()
        {
            IFetiDP3dMatrixManagerFactory matricesFactory = new FetiDP3dMatrixManagerFactoryDense();
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            FetiDP3dMatrixManagerSerial matrixManager = PrepareCoarseProblemSubdomainMatrices(model, dofSeparator,
                lagrangesEnumerator, augmentationConstraints, matricesFactory, MatrixFormat.Dense);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Input data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.LinearSystem.RhsConcrete = ExpectedSubdomainMatrices.GetVectorFf(sub.ID);

                // Prepare the subdomain data
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
            }

            // Calculate the global data to test
            matrixManager.CalcCoarseProblemRhs();

            // Check
            double tol = 1E-13;
            Assert.True(ExpectedGlobalMatrices.VectorDr.Equals(matrixManager.GlobalDr, tol));
        }

        [Fact]
        public static void TestVectorsFbcFr()
        {
            IFetiDP3dMatrixManagerFactory matricesFactory = new FetiDP3dMatrixManagerFactoryDense();
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            FetiDP3dMatrixManagerSerial matrixManager = new FetiDP3dMatrixManagerSerial(model, dofSeparator, lagrangesEnumerator,
                augmentationConstraints, matricesFactory);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Calculate the necessary vectors
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.LinearSystem.RhsConcrete = ExpectedSubdomainMatrices.GetVectorFf(sub.ID);
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();

                // Check
                double tol = 1E-13;
                Assert.True(ExpectedSubdomainMatrices.GetVectorFbc(sub.ID).Equals(subdomainMatrices.Fbc, tol));
                Assert.True(ExpectedSubdomainMatrices.GetVectorFr(sub.ID).Equals(subdomainMatrices.Fr, tol));
            }
        }

        internal static FetiDP3dMatrixManagerSerial PrepareCoarseProblemSubdomainMatrices(IModel model,
            IFetiDPDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangesEnumerator, 
            IAugmentationConstraints augmentationConstraints, IFetiDP3dMatrixManagerFactory matricesFactory, MatrixFormat format)
        {
            var matrixManager = new FetiDP3dMatrixManagerSerial(model, dofSeparator, lagrangesEnumerator, 
                augmentationConstraints, matricesFactory);
            SetKffRhs(model, matrixManager, format);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Prepare the subdomain data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
                subdomainMatrices.InvertKrr(true);
            }

            return matrixManager;
        }

        internal static void PrepareCoarseProblemSubdomainMatrices(IModel model, FetiDP3dMatrixManagerSerial matrixManager, 
            MatrixFormat format)
        {
            SetKffRhs(model, matrixManager, format);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Prepare the subdomain data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                subdomainMatrices.ExtractCornerRemainderSubmatrices();
                subdomainMatrices.ExtractCornerRemainderRhsSubvectors();
                subdomainMatrices.InvertKrr(true);
            }
        }

        private static void SetKffRhs(IModel model, IFetiDPMatrixManager matrixManager, MatrixFormat format)
        {
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                // Input data
                IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(sub);
                Matrix Kff = ExpectedSubdomainMatrices.GetMatrixKff(sub.ID);
                if ((format == MatrixFormat.Dense) || (format == MatrixFormat.Skyline))
                {
                    var castedLS = (SingleSubdomainSystemMpi<SkylineMatrix>)subdomainMatrices.LinearSystem;
                    castedLS.Matrix = SkylineMatrix.CreateFromMatrix(Kff);
                }
                else
                {
                    var castedLS = (SingleSubdomainSystemMpi<DokSymmetric>)subdomainMatrices.LinearSystem;
                    castedLS.Matrix = DokSymmetric.CreateFromMatrix(Kff);
                }
                subdomainMatrices.LinearSystem.RhsConcrete = ExpectedSubdomainMatrices.GetVectorFf(sub.ID);
            }
        }
    }
}
