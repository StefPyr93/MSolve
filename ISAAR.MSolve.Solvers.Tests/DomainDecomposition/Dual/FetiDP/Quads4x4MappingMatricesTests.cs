using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4MappingMatricesTests
    {
        [Fact] // TODO: Can be removed as it tests old code and the replacement is tested elsewhere
        public static void TestDofSeparation()
        {
            // Create model
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparatorOLD();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> cornerNodes = Example4x4QuadsHomogeneous.DefineCornerNodesSubdomain(subdomain);
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, cornerNodes);
            }
            
            // Check
            for (int s = 0; s < 4; ++s)
            {
                (int[] cornerDofs, int[] remainderDofs, int[] boundaryRemainderDofs, int[] internalRemainderDofs) =
                    Example4x4QuadsHomogeneous.GetDofSeparation(s);
                ArrayChecking.CheckEqual(cornerDofs, dofSeparator.CornerDofIndices[s]);
                ArrayChecking.CheckEqual(remainderDofs, dofSeparator.RemainderDofIndices[s]);
                ArrayChecking.CheckEqual(boundaryRemainderDofs, dofSeparator.BoundaryDofIndices[s]);
                ArrayChecking.CheckEqual(internalRemainderDofs, dofSeparator.InternalDofIndices[s]);
            }
        }

        [Fact] // TODO: Can be removed as it tests old code and the replacement is tested elsewhere
        public static void TestSignedBooleanMatrices()
        {
            // Create model
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparatorOLD();
            dofSeparator.DefineGlobalBoundaryDofs(model, Example4x4QuadsHomogeneous.DefineCornerNodesGlobal(model));
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> cornerNodes = Example4x4QuadsHomogeneous.DefineCornerNodesSubdomain(subdomain);
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, cornerNodes);
            }

            // Enumerate lagranges
            var crosspointStrategy = new FullyRedundantConstraints();
            var lagrangeEnumerator = new FetiDPLagrangeMultipliersEnumeratorOLD(crosspointStrategy, dofSeparator);
            lagrangeEnumerator.DefineBooleanMatrices(model);

            // Check
            int expectedNumLagrangeMultipliers = 8;
            Assert.Equal(expectedNumLagrangeMultipliers, lagrangeEnumerator.NumLagrangeMultipliers);
            double tolerance = 1E-13;
            for (int s = 0; s < 4; ++s)
            {
                Matrix Br = lagrangeEnumerator.BooleanMatrices[s].CopyToFullMatrix(false);
                Matrix expectedBr = Example4x4QuadsHomogeneous.GetMatrixBr(s);
                Assert.True(expectedBr.Equals(Br, tolerance));
            }
        }

        [Fact] // TODO: Can be removed as it tests old code and the replacement is tested elsewhere
        public static void TestUnsignedBooleanMatrices()
        {
            // Create model
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparatorOLD();
            dofSeparator.DefineGlobalCornerDofs(model, Example4x4QuadsHomogeneous.DefineCornerNodesGlobal(model));
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> cornerNodes = Example4x4QuadsHomogeneous.DefineCornerNodesSubdomain(subdomain);
                //IEnumerable<INode> remainderAndConstrainedNodes = subdomain.Nodes.Where(node => !cornerNodes[s].Contains(node));
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, cornerNodes);
            }
            dofSeparator.CalcCornerMappingMatrices(model);

            // Check
            int expectedNumCornerDofs = 8;
            Assert.Equal(expectedNumCornerDofs, dofSeparator.NumGlobalCornerDofs);
            double tolerance = 1E-13;
            for (int s = 0; s < 4; ++s)
            {
                UnsignedBooleanMatrix Bc = dofSeparator.CornerBooleanMatrices[s];
                Matrix expectedBc = Example4x4QuadsHomogeneous.GetMatrixBc(s);
                Assert.True(expectedBc.Equals(Bc, tolerance));
            }
        }
    }
}
