using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.XFEM.FEM.Mesh
{
    public class Vertex
    {
        public Vertex(int id, params double[] coords)
        {
            ID = id;
            Coords = coords;
        }

        public int ID { get; set; }

        public double[] Coords { get; set; }
    }
}
