using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using MPI;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests.Transfer
{
    internal static class SubdomainComparisons
    {
        internal static void CheckSameElements(XSubdomain expectedSubdomain, XSubdomain actualSubdomain)
        {
            Assert.Equal(expectedSubdomain.NumElements, actualSubdomain.NumElements);
            foreach (IXFiniteElement expectedElement in expectedSubdomain.Elements.Values)
            {
                Assert.True(actualSubdomain.Elements.ContainsKey(expectedElement.ID));
                IXFiniteElement actualElement = actualSubdomain.Elements[expectedElement.ID];
                Assert.Equal(expectedElement.ID, actualElement.ID);
                Assert.Equal(expectedElement.Nodes.Count, actualElement.Nodes.Count);
                for (int n = 0; n < expectedElement.Nodes.Count; ++n)
                {
                    Assert.Equal(expectedElement.Nodes[n].ID, actualElement.Nodes[n].ID);
                }
            }
        }

        internal static void CheckSameNodalDisplacements(XSubdomain expectedSubdomain, XSubdomain actualSubdomain)
        {
            int precision = 8;
            foreach (XNode expectedNode in expectedSubdomain.Nodes.Values)
            {
                List<Constraint> expectedDisplacements = expectedNode.Constraints;
                List<Constraint> actualDisplacements = actualSubdomain.Nodes[expectedNode.ID].Constraints;
                Assert.Equal(expectedDisplacements.Count, actualDisplacements.Count);
                for (int d = 0; d < expectedDisplacements.Count; ++d)
                {
                    Assert.Equal(expectedDisplacements[d].DOF, actualDisplacements[d].DOF);
                    Assert.Equal(expectedDisplacements[d].Amount, actualDisplacements[d].Amount, precision);
                }
            }
        }

        internal static void CheckSameNodalLoads(XSubdomain expectedSubdomain, XSubdomain actualSubdomain)
        {
            int precision = 8;
            Assert.Equal(expectedSubdomain.NumNodalLoads, actualSubdomain.NumNodalLoads);
            for (int nl = 0; nl < expectedSubdomain.NumNodalLoads; ++nl)
            {
                NodalLoad expectedLoad = expectedSubdomain.NodalLoads[nl];
                NodalLoad actualLoad = actualSubdomain.NodalLoads[nl];
                Assert.Equal(expectedLoad.Node.ID, actualLoad.Node.ID);
                Assert.Equal(expectedLoad.DofType, actualLoad.DofType);
                Assert.Equal(expectedLoad.Amount, actualLoad.Amount, precision);
            }
        }

        internal static void CheckSameNodes(XSubdomain expectedSubdomain, XSubdomain actualSubdomain)
        {
            int precision = 8;
            Assert.Equal(expectedSubdomain.NumNodes, actualSubdomain.NumNodes);
            foreach (XNode expectedNode in expectedSubdomain.Nodes.Values)
            {
                Assert.True(actualSubdomain.Nodes.ContainsKey(expectedNode.ID));
                XNode actualNode = actualSubdomain.Nodes[expectedNode.ID];
                Assert.Equal(expectedNode.ID, actualNode.ID);
                Assert.Equal(expectedNode.X, actualNode.X, precision);
                Assert.Equal(expectedNode.Y, actualNode.Y, precision);
                Assert.Equal(expectedNode.Z, actualNode.Z, precision);
            }
        }
    }
}
