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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests
{
    public static class FetiDPDisplacementsCalculatorSerialTests
    {
        [Fact]
        public static void TestCornerDisplacements()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(Example4x4QuadsHomogeneous.CreateModel());
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            Vector lagranges = Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers;
            var displacementsCalculator = new FreeDofDisplacementsCalculatorSerial(model, dofSeparator, matrixManager,
                lagrangesEnumerator);

            MethodInfo method = displacementsCalculator.GetType().GetMethod("CalcCornerDisplacements",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector uc = (Vector)method.Invoke(displacementsCalculator, new object[] { flexibility, lagranges });

            double tol = 1E-12;
            Assert.True(Example4x4QuadsHomogeneous.SolutionCornerDisplacements.Equals(uc, tol));
        }

        [Fact]
        public static void TestFreeDisplacements()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();

            var displacementsCalculator = new FreeDofDisplacementsCalculatorSerial(model, dofSeparator, matrixManager, 
                lagrangesEnumerator);
            Vector lagranges = Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers;
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);

            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                double tol = 1E-7;
                IVectorView uf = matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.Solution;
                Assert.True(Example4x4QuadsHomogeneous.GetSolutionFreeDisplacements(sub.ID).Equals(uf, tol));
            }
            
        }
    }
}
