using System;
using System.Collections.Generic;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dDisplacementsCalculatorSerialTests
    {
        [Fact]
        public static void TestCornerDisplacementsAndAugmentedLagranges()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            Vector lagranges = Example4x4x4Quads.ExpectedSolutions.SolutionLagrangesSimple();

            var displacementsCalculator = new FetiDP3dFreeDofDisplacementsCalculatorSerial(model, dofSeparator, 
                lagrangesEnumerator, augmentationConstraints, matrixManager); 

            MethodInfo method = displacementsCalculator.GetType().GetMethod("CalcCornerDisplacementsAndAugmentedLagranges",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector ucTilde = (Vector)method.Invoke(displacementsCalculator, new object[] { flexibility, lagranges });

            double tol = 1E-4;
            Assert.True(Example4x4x4Quads.ExpectedSolutions.ucTilde().Equals(ucTilde, tol));

        }

        [Fact]
        public static void TestFreeDisplacements()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, dofSeparator, lagrangesEnumerator);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            Vector lagranges = Example4x4x4Quads.ExpectedSolutions.SolutionLagrangesSimple();

            var displacementsCalculator = new FetiDP3dFreeDofDisplacementsCalculatorSerial(model, dofSeparator, lagrangesEnumerator, augmentationConstraints, matrixManager);
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);
           

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                double tol = 1E-4;
                IVectorView uf = matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.Solution;
                Assert.True(Example4x4x4Quads.ExpectedSolutions.uFree(sub.ID).Equals(uf, tol));
            }
            
        }
    }
}
