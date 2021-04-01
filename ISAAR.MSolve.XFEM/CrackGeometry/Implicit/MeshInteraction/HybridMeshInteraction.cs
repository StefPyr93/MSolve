using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    class HybridMeshInteraction : IMeshInteraction
    {

        public HybridMeshInteraction()
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

            //Warning: this might actually be worse than Stolarska's criterion. At least that one enriched the dubious 
            //intersected elements with tip enrichments. This one just ignores them.
            if (minBodyLevelSet * maxBodyLevelSet <= 0.0)
            {
                if (minTipLevelSet * maxTipLevelSet <= 0)
                {
                    var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                    if (outline.IsPointInsidePolygon(crackTip)) return CrackElementPosition.ContainsTip;
                }
                else if (maxTipLevelSet < 0) return CrackElementPosition.Intersected;
            }
            return CrackElementPosition.Irrelevant;
        }
    }
}
