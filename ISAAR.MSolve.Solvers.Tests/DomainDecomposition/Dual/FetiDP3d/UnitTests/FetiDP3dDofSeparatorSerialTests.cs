using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dDofSeparatorSerialTests
    {
        [Fact]
        public static void TestDofSeparation()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator) = CreateModelAndDofSeparator();

            // Check
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {

                int[] cornerDofs = ExpectedConnectivityData.GetDofsCornerInFree(sub.ID);
                int[] remainderDofs = ExpectedConnectivityData.GetDofsRemainderInFree(sub.ID);
                int[] boundaryRemainderDofs = ExpectedConnectivityData.GetDofsBoundaryInRemainder(sub.ID);
                int[] internalRemainderDofs = ExpectedConnectivityData.GetDofsInternalInRemainder(sub.ID);
                    
                ArrayChecking.CheckEqual(cornerDofs, dofSeparator.GetCornerDofIndices(sub));
                ArrayChecking.CheckEqual(remainderDofs, dofSeparator.GetRemainderDofIndices(sub));
                ArrayChecking.CheckEqual(boundaryRemainderDofs, dofSeparator.GetBoundaryDofIndices(sub));
                ArrayChecking.CheckEqual(internalRemainderDofs, dofSeparator.GetInternalDofIndices(sub));
            }
        }

        [Fact]
        public static void TestCornerBooleanMatrices()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator) = CreateModelAndDofSeparator();

            // Check
            int expectedNumCornerDofs = 3;
            Assert.Equal(expectedNumCornerDofs, dofSeparator.NumGlobalCornerDofs);
            double tolerance = 1E-13;
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(sub);
                Matrix expectedBc = ExpectedConnectivityData.GetMatrixBc(sub.ID);
                Assert.True(expectedBc.Equals(Bc, tolerance));
            }
        }

        internal static (IModel, FetiDPDofSeparatorSerial) CreateModelAndDofSeparator()
        {
            // Create model
            Model model = ModelCreator.CreateModel();
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Separate dofs and calculate the boolean matrices
            ICornerNodeSelection cornerNodes = ModelCreator.DefineCornerNodeSelectionSerial(model);
            var dofSeparator = new FetiDPDofSeparatorSerial(model, cornerNodes);
            dofSeparator.SeparateDofs(new MockSeparatedDofReordering());

            return (model, dofSeparator);
        }
    }
}
