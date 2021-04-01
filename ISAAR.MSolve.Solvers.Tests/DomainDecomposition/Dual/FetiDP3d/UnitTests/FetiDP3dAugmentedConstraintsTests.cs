using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;
using Xunit;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dAugmentedConstraintsTests
    {
        private enum ConstraintsSelection
        {
            Simple, Redundant 
        }

        [Theory]
        [InlineData(ConstraintsSelection.Simple)]
        [InlineData(ConstraintsSelection.Redundant)]
        private static void TestMatrixQr(ConstraintsSelection constraints)
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            IAugmentationConstraints augmentationConstraints;
            int numAugmentationConstraints; 
            Matrix QrExpected;
            if (constraints == ConstraintsSelection.Simple)
            {
                augmentationConstraints = CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
                numAugmentationConstraints = ExpectedConnectivityData.NumGlobalAugmentationConstraintsSimple;
                QrExpected = ExpectedConnectivityData.MatrixQrSimple;
            }
            else
            {
                augmentationConstraints = CalcAugmentationConstraintsRedundant(model, lagrangesEnumerator);
                numAugmentationConstraints = ExpectedConnectivityData.NumGlobalAugmentationConstraintsRedundant;
                QrExpected = ExpectedConnectivityData.MatrixQrRedundant;
            }

            //string pathEx = @"C:\Users\Serafeim\Desktop\MPI\Tests\expectedQr.txt";
            //string path = @"C:\Users\Serafeim\Desktop\MPI\Tests\Qr.txt";
            //string pathDiff = @"C:\Users\Serafeim\Desktop\MPI\Tests\diffQr.txt";
            //var writer = new LinearAlgebra.Output.FullMatrixWriter();
            //writer.WriteToFile(QrExpected, pathEx);
            //writer.WriteToFile(augmentationConstraints.MatrixQr, path);
            //writer.WriteToFile(QrExpected - augmentationConstraints.MatrixQr, pathDiff);

            // Check
            Assert.Equal(numAugmentationConstraints, augmentationConstraints.NumGlobalAugmentationConstraints);
            double tolerance = 1E-13;
            IMappingMatrix Qr = augmentationConstraints.MatrixGlobalQr;
            Matrix explicitQr = Qr.MultiplyRight(Matrix.CreateIdentity(Qr.NumColumns));
            Assert.True(QrExpected.Equals(explicitQr, tolerance));
        }

        [Fact]
        private static void TestMatricesQ1Ba()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();

            IAugmentationConstraints augmentationConstraints = 
                CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);

            //Matrix expectedQ = augmentationConstraints.MatrixGlobalQr;
            Matrix expectedQ = ExpectedConnectivityData.MatrixQrSimple;

            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                Matrix R1 = augmentationConstraints.GetMatrixR1(subdomain).CopyToFullMatrix();
                IMappingMatrix Ba = augmentationConstraints.GetMatrixBa(subdomain);
                Matrix explicitBa = Ba.MultiplyRight(Matrix.CreateIdentity(Ba.NumColumns));

                Matrix computedBrTransposeTimesQr = R1 * Matrix.CreateFromMatrix(explicitBa);

                Matrix expectedBr = ExpectedConnectivityData.GetMatrixBr(subdomain.ID);
                Matrix expectedBrTransposeTimesQr = expectedBr.MultiplyRight(expectedQ, true);

                // Check
                double tolerance = 1E-13;
                Assert.True(expectedBrTransposeTimesQr.Equals(computedBrTransposeTimesQr, tolerance));
            }
        }

        internal static IAugmentationConstraints CalcAugmentationConstraintsSimple(IModel model, IFetiDPDofSeparator dofSeparator,
            LagrangeMultipliersEnumeratorSerial lagrangesEnumerator)
        {
            Dictionary<ISubdomain, HashSet<INode>> midsideNodes = ModelCreator.DefineMidsideNodesAll(model);
            IMidsideNodesSelection midsideNodesSelection = new UserDefinedMidsideNodes(midsideNodes,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
            
            IAugmentationConstraints augmentationConstraints =
                    new AugmentationConstraints(model, midsideNodesSelection, dofSeparator, lagrangesEnumerator);
            augmentationConstraints.CalcAugmentationMappingMatrices();
            return augmentationConstraints;
        }

        internal static IAugmentationConstraints CalcAugmentationConstraintsRedundant(IModel model,
            LagrangeMultipliersEnumeratorSerial lagrangesEnumerator)
        {
            Dictionary<ISubdomain, HashSet<INode>> midsideNodes = ModelCreator.DefineMidsideNodesAll(model);
            IMidsideNodesSelection midsideNodesSelection = new UserDefinedMidsideNodes(midsideNodes,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
            IAugmentationConstraints augmentationConstraints =
                    new AugmentationConstraintsRedundant(model, midsideNodesSelection, lagrangesEnumerator);
            augmentationConstraints.CalcAugmentationMappingMatrices();
            return augmentationConstraints;
        }
    }
}
