using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    interface IMeshInteraction
    {
        CrackElementPosition FindRelativePositionOf(XContinuumElement2D element, CartesianPoint crackTip,
            Dictionary<XNode, double> levelSetsBody, Dictionary<XNode, double> levelSetsTip);
    }
}
