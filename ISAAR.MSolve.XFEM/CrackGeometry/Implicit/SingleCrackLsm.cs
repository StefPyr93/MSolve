using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit
{
    public class SingleCrackLsm
    {
        private static readonly IComparer<CartesianPoint> pointComparer = new Point2DComparerXMajor();

        private readonly Dictionary<XNode, double> levelSetsBody;
        private readonly Dictionary<XNode, double> levelSetsTip;
        private readonly IMeshInteraction meshInteraction;
        

        public SingleCrackLsm(Dictionary<XNode, double> levelSetsBody, Dictionary<XNode, double> levelSetsTip,
            CartesianPoint crackTip)
        {
            this.CrackTip = crackTip;
            this.levelSetsBody = levelSetsBody;
            this.levelSetsTip = levelSetsTip;
            this.meshInteraction = new SerafeimMeshInteraction();
        }
        public CartesianPoint CrackTip { get; }

        public Dictionary<XNode, double> LevelSetsBody => levelSetsBody;
        public Dictionary<XNode, double> LevelSetsTip => levelSetsTip;

        public SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element)
        {
            var triangleVertices = new SortedSet<CartesianPoint>(element.Nodes, pointComparer);
            int nodesCount = element.Nodes.Count;
            CrackElementPosition relativePosition = meshInteraction.FindRelativePositionOf(element, CrackTip,
                levelSetsBody, levelSetsTip);
            if (relativePosition != CrackElementPosition.Irrelevant)
            {
                // Find the intersections between element edges and the crack. TODO: See Serafeim's Msc Thesis for a correct procedure.
                for (int i = 0; i < nodesCount; ++i)
                {
                    XNode node1 = element.Nodes[i];
                    XNode node2 = element.Nodes[(i + 1) % nodesCount];
                    double levelSet1 = levelSetsBody[node1];
                    double levelSet2 = levelSetsBody[node2];

                    if (levelSet1 * levelSet2 < 0.0)
                    {
                        // The intersection point between these nodes can be found using the linear interpolation, see 
                        // Sukumar 2001
                        double k = -levelSet1 / (levelSet2 - levelSet1);
                        double x = node1.X + k * (node2.X - node1.X);
                        double y = node1.Y + k * (node2.Y - node1.Y);

                        // TODO: For the tip element one intersection point is on the crack extension and does not  
                        // need to be added. It is not wrong though.
                        triangleVertices.Add(new CartesianPoint(x, y));
                    }
                    else if (levelSet1 == 0.0) triangleVertices.Add(node1); // TODO: perhaps some tolerance is needed.
                    else if (levelSet2 == 0.0) triangleVertices.Add(node2);
                }

                if (relativePosition == CrackElementPosition.ContainsTip) triangleVertices.Add(CrackTip);
            }
            return triangleVertices;
        }

        // TODO: with narrow band this should throw an exception if the node is not tracked.
        public double SignedDistanceOf(XNode node) => levelSetsBody[node];

        // TODO: with narrow band this should throw an exception if the element is not tracked.
        public double SignedDistanceOf(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation)
        {
            double signedDistance = 0.0;
            for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
            {
                signedDistance += interpolation.ShapeFunctions[nodeIdx] * levelSetsBody[element.Nodes[nodeIdx]];
            }
            return signedDistance;
        }

        public Tuple<double, double> SignedDistanceGradientThrough(NaturalPoint point,
            IReadOnlyList<XNode> elementNodes, EvalInterpolation2D interpolation)
        {
            double gradientX = 0.0;
            double gradientY = 0.0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                double dNdx = interpolation.ShapeGradientsCartesian[nodeIdx, 0];
                double dNdy = interpolation.ShapeGradientsCartesian[nodeIdx, 1];

                double levelSet = levelSetsBody[elementNodes[nodeIdx]];
                gradientX += dNdx * levelSet;
                gradientY += dNdy * levelSet;
            }
            return new Tuple<double, double>(gradientX, gradientY);
        }
    }
}
