using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Fill in connectivity mappings for wedges and pyramids. The documentation for gmsh 3.0 has nice figures.
namespace MGroup.XFEM.FEM.Mesh.GMSH
{
    /// <summary>
    /// Converts cell types and the order of their vertices from GMSH to MSolve. WARNING: These cell codes are for files 
    /// following the 2.2 mesh format of GMSH. Other formats may differ.
    /// </summary>
    internal class GmshCellFactory
    {
        /// <summary>
        /// Keys are the GMSH cell types. Values are the MSolve mesh types.
        /// </summary>
        private static readonly IReadOnlyDictionary<int, CellType> cellTypes;

        /// <summary>
        /// Vertex order for cells. Index = gmsh order, value = MSolve order.
        /// </summary>
        private static readonly IReadOnlyDictionary<CellType, int[]> gmshToMSolveVertexOrders;

        static GmshCellFactory()
        {
            var codes = new Dictionary<int, CellType>();
            codes[1] = CellType.Line;
            codes[2] = CellType.Tri3;
            codes[3] = CellType.Quad4;
            codes[4] = CellType.Tet4;
            codes[5] = CellType.Hexa8;
            codes[6] = CellType.Wedge6;
            codes[7] = CellType.Pyra5;
            //codes[8] = CellType.Line3;
            codes[9] = CellType.Tri6;
            codes[10] = CellType.Quad9;
            codes[11] = CellType.Tet10;
            codes[12] = CellType.Hexa27;
            codes[13] = CellType.Wedge18;
            codes[14] = CellType.Pyra14;
            //codes[15] = CellType.SinglePoint;
            codes[16] = CellType.Quad8;
            codes[17] = CellType.Hexa20;
            codes[18] = CellType.Wedge15;
            codes[19] = CellType.Pyra13;
            cellTypes = codes;

            var connectivity = new Dictionary<CellType, int[]>();
            connectivity[CellType.Line] = new int[] { 0, 1 };
            connectivity[CellType.Tri3] = new int[] { 2, 0, 1 };
            connectivity[CellType.Quad4] = new int[] { 0, 1, 2, 3 };
            connectivity[CellType.Tet4] = new int[] { 0, 1, 2, 3 };
            connectivity[CellType.Hexa8] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            connectivity[CellType.Wedge6] = null;
            connectivity[CellType.Pyra5] = null;
            //connectivity[CellType.Line3] = new int[] { 0, 1, 2 };
            connectivity[CellType.Tri6] = new int[] { 2, 0, 1, 5, 3, 4 };
            connectivity[CellType.Quad9] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
            connectivity[CellType.Tet10] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            connectivity[CellType.Hexa27] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 };
            connectivity[CellType.Wedge18] = null;
            connectivity[CellType.Pyra14] = null;
            //connectivity[CellType.SinglePoint] = new int[] { 0 };
            connectivity[CellType.Quad8] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            connectivity[CellType.Hexa20] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
            connectivity[CellType.Wedge15] = null;
            connectivity[CellType.Pyra13] = null;
            gmshToMSolveVertexOrders = connectivity;
        }

        private readonly IReadOnlyList<Vertex> allVertices;

        public GmshCellFactory(IReadOnlyList<Vertex> allVertices)
        {
            this.allVertices = allVertices;
        }

        public HashSet<int> GetValidGmshCellTypes(int minDimension)
        {
            int[] cellTypes1D = { 1, 8 };
            int[] cellTypes2D = { 2, 3, 9, 10, 16 };
            int[] cellTypes3D = { 4, 5, 6, 7, 11, 12, 13, 14, 17, 18, 19 };

            var validCellTypes = new HashSet<int>(cellTypes3D);
            if (minDimension <= 2) validCellTypes.UnionWith(cellTypes2D);
            if (minDimension == 1) validCellTypes.UnionWith(cellTypes1D);

            return validCellTypes;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="cellID"></param>
        /// <param name="gmshCellType"></param>
        /// <param name="gmshVertexIDs">Must be 0-based.</param>
        /// <param name="physicalGroup"></param>
        /// <returns></returns>
        public Cell CreateGmshCell(int cellID, int gmshCellType, int[] gmshVertexIDs, int physicalGroup)
        {
            CellType msolveCellType = cellTypes[gmshCellType];
            int[] gmshToMSolveVertices = gmshToMSolveVertexOrders[msolveCellType];

            // GMSH to MSolve order of vertices
            int numVertices = gmshVertexIDs.Length;
            Debug.Assert(gmshToMSolveVertices.Length == numVertices);
            var msolveVertexIDs = new int[numVertices];
            for (int i = 0; i < numVertices; ++i)
            {
                msolveVertexIDs[gmshToMSolveVertices[i]] = gmshVertexIDs[i];
            }

            // Possibly fix the order of vertices to avoid negative Jacobians
            //TODO: Needs something for Hexa8 at least. For 2nd order elements, I cannot just approximate with polygons at their
            //      nodes, since these polygons will not always be convex.
            if ((msolveCellType == CellType.Tri3) || (msolveCellType == CellType.Quad4))
            {
                FixVertexOrderTri3Quad4(msolveVertexIDs);
            }
            else if (msolveCellType == CellType.Tet4)
            {
                FixVertexOrderTet4(msolveVertexIDs);
            }

            return new Cell() 
            { 
                ID = cellID, CellType = msolveCellType, VertexIDs = msolveVertexIDs, PhysicalGroup = physicalGroup 
            };
        }

        private void FixVertexOrderTri3Quad4(int[] vertexIDs)
        {
            List<double[]> vertices = SelectVertices(vertexIDs);

            // Calculate the sign of the area (here double the area)
            double cellArea = 0.0; 
            for (int i = 0; i < vertices.Count; ++i)
            {
                double[] vertex0 = vertices[i];
                double[] vertex1 = vertices[(i + 1) % vertices.Count];
                cellArea += vertex0[0] * vertex1[1] - vertex1[0] * vertex0[1];
            }
            int signOfArea = Math.Sign(cellArea);
            
            if (signOfArea < 0) Array.Reverse(vertexIDs); // The area of the cell with clockwise vertices is negative!
            else if (signOfArea == 0) throw new Exception("Degenerate cell.");
        }

        private void FixVertexOrderTet4(int[] vertexIDs)
        {
            List<double[]> vertices = SelectVertices(vertexIDs);

            // Vectors from origin to vertices
            var p0 = Vector.CreateFromArray(vertices[0]);
            var p1 = Vector.CreateFromArray(vertices[1]);
            var p2 = Vector.CreateFromArray(vertices[2]);
            var p3 = Vector.CreateFromArray(vertices[3]);

            // Vectors along the edges of the Tet4. The order is important.
            Vector v01 = p1 - p0;
            Vector v02 = p2 - p0;
            Vector v03 = p3 - p0;

            // If the normal vector of the face p0-p1-p2 points towards p3, then the vertex order is ok.
            Vector n = v01.CrossProduct(v02);
            int sign = Math.Sign(v03 * n);
            if (sign < 0)
            {
                // The normal vector of the face p0-p1-p2 points away from p3, but the normal of the face p2-p1-p0 points 
                // towards p3. Swap p0 and p2.
                int swap = vertexIDs[0];
                vertexIDs[0] = vertexIDs[2];
                vertexIDs[2] = swap;
            }
            else if (sign == 0) throw new Exception("Degenerate cell.");
        }

        private List<double[]> SelectVertices(int[] vertexIDs)
        {
            var vertices = new List<double[]>(vertexIDs.Length);
            for (int i = 0; i < vertexIDs.Length; ++i)
            {
                //TODO: This only works if they are not contiguous. GmshReader could enforce this.
                vertices.Add(allVertices[vertexIDs[i]].Coords);
            }
            return vertices;
        }
    }
}
