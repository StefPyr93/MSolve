using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Mesh;

namespace MGroup.XFEM.FEM.Mesh
{
    public class Cell
    {
        public int ID { get; set; }

        public CellType CellType { get; set; }

        public int[] VertexIDs { get; set; }

        public int PhysicalGroup { get; set; } = 0;
    }
}
