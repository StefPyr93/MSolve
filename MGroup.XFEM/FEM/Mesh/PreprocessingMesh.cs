using System;
using System.Collections.Generic;
using System.Text;

//TODO: perhaps this should not be just for preprocessing purposes. Instead it could have mesh operations (e.g. find neighboring 
//      elements), partitioning operations and any connectivity feature that does not depend on FEM, XFEM, IGA, materials, etc.
namespace MGroup.XFEM.FEM.Mesh
{
    public class PreprocessingMesh
    {
        public List<Cell> Cells { get; set; }

        public List<Vertex> Vertices { get; set; }
    }
}
