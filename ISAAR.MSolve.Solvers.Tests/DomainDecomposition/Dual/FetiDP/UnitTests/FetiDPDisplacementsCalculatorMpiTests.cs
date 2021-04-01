using System;
using System.Collections.Generic;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
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
    public static class FetiDPDisplacementsCalculatorMpiTests
    {
        public static void TestCornerDisplacements(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator,
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(Example4x4QuadsHomogeneous.CreateModel());
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();
            Vector lagranges = Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers;
            var displacementsCalculator = new FreeDofDisplacementsCalculatorMpi(procs, model, dofSeparator, matrixManager,
                lagrangesEnumerator);

            MethodInfo method = displacementsCalculator.GetType().GetMethod("CalcCornerDisplacements",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector uc = (Vector)method.Invoke(displacementsCalculator, new object[] { flexibility, lagranges });

            if (procs.IsMasterProcess)
            {
                double tol = 1E-12;
                Assert.True(Example4x4QuadsHomogeneous.SolutionCornerDisplacements.Equals(uc, tol));
            }
        }

        public static void TestFreeDisplacements(int numProcesses)
        {
            (ProcessDistribution procs, IModel model, FetiDPDofSeparatorMpi dofSeparator, 
                LagrangeMultipliersEnumeratorMpi lagrangesEnumerator) =
                FetiDPLagrangesEnumeratorMpiTests.CreateModelDofSeparatorLagrangesEnumerator(numProcesses);
            IFetiDPMatrixManager matrixManager = new MockMatrixManager(model);
            IFetiDPFlexibilityMatrix flexibility = new MockFlexibilityMatrix();

            var displacementsCalculator = new FreeDofDisplacementsCalculatorMpi(procs, model, dofSeparator, matrixManager, 
                lagrangesEnumerator);
            Vector lagranges = Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers;
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);


            double tol = 1E-7;
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                IVectorView uf = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).LinearSystem.Solution;
                Assert.True(Example4x4QuadsHomogeneous.GetSolutionFreeDisplacements(subdomain.ID).Equals(uf, tol));
            }
        }
    }
}
