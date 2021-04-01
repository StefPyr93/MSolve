using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    class StolarskaMeshInteraction: IMeshInteraction
    {
        public StolarskaMeshInteraction()
        {
        }

        public CrackElementPosition FindRelativePositionOf(XContinuumElement2D element, CartesianPoint crackTip,
            Dictionary<XNode, double> levelSetsBody, Dictionary<XNode, double> levelSetsTip)
        {
            double minBodyLevelSet = double.MaxValue;
            double maxBodyLevelSet = double.MinValue;
            double minTipLevelSet = double.MaxValue;
            double maxTipLevelSet = double.MinValue;

            foreach (XNode node in element.Nodes)
            {
                double bodyLevelSet = levelSetsBody[node];
                double tipLevelSet = levelSetsTip[node];
                if (bodyLevelSet < minBodyLevelSet) minBodyLevelSet = bodyLevelSet;
                if (bodyLevelSet > maxBodyLevelSet) maxBodyLevelSet = bodyLevelSet;
                if (tipLevelSet < minTipLevelSet) minTipLevelSet = tipLevelSet;
                if (tipLevelSet > maxTipLevelSet) maxTipLevelSet = tipLevelSet;
            }

            // Warning: This criterion might give false positives for tip elements (see Serafeim's thesis for details)
            if (minBodyLevelSet * maxBodyLevelSet <= 0.0)
            {
                if (minTipLevelSet * maxTipLevelSet <= 0) return CrackElementPosition.ContainsTip;
                else if (maxTipLevelSet < 0) return CrackElementPosition.Intersected;
            }
            return CrackElementPosition.Irrelevant;
        }
    }
}
