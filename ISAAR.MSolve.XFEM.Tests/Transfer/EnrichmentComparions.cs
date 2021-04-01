using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests.Transfer
{
    internal static class EnrichmentComparions
    {
        internal static void CheckSameEnrichments(XSubdomain expectedSubdomain, XSubdomain actualSubdomain, int precision = 8)
        {
            //TODO: Unfortunately I cannot compare the enrichments based on IDs. Perhaps they should have

            // Node enrichments and their nodal values.
            foreach (XNode expectedNode in expectedSubdomain.Nodes.Values)
            {
                var expectedEnrichments = new List<double>();
                foreach (double[] enrichments in expectedNode.EnrichmentItems.Values)
                {
                    foreach (double val in enrichments) expectedEnrichments.Add(val);
                }

                XNode actualNode = actualSubdomain.Nodes[expectedNode.ID];
                var actualEnrichments = new List<double>();
                foreach (double[] enrichments in actualNode.EnrichmentItems.Values)
                {
                    foreach (double val in enrichments) actualEnrichments.Add(val);
                }

                Assert.Equal(expectedNode.EnrichmentItems.Count, actualNode.EnrichmentItems.Count);
                Assert.Equal(expectedEnrichments.Count, actualEnrichments.Count);
                for (int i = 0; i < expectedEnrichments.Count; ++i)
                {
                    Assert.Equal(expectedEnrichments[i], actualEnrichments[i], precision);
                }
            }

            // Element enrichments
            foreach (IXFiniteElement expectedElement in expectedSubdomain.Elements.Values)
            {
                IXFiniteElement actualElement = actualSubdomain.Elements[expectedElement.ID];
                Assert.Equal(expectedElement.EnrichmentItems.Count, actualElement.EnrichmentItems.Count);
            }
        }

        internal static void CheckSameLevelSets(SingleCrackLsm expectedLsm, XSubdomain expectedSubdomain,
            SingleCrackLsm actualLsm, XSubdomain actualSubdomain, int precision = 8)
        {
            foreach (XNode expectedNode in expectedSubdomain.Nodes.Values)
            {
                double expectedBodyLevelSet = expectedLsm.LevelSetsBody[expectedNode];
                double expectedTipLevelSet = expectedLsm.LevelSetsBody[expectedNode];

                XNode actualNode = actualSubdomain.Nodes[expectedNode.ID];
                double actualBodyLevelSet = actualLsm.LevelSetsBody[actualNode];
                double actualTipLevelSet = actualLsm.LevelSetsBody[actualNode];

                Assert.Equal(expectedBodyLevelSet, actualBodyLevelSet, precision);
                Assert.Equal(expectedTipLevelSet, actualTipLevelSet, precision);
            }
        }
    }
}
