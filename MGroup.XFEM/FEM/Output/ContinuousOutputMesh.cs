using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using MGroup.XFEM.Output.Mesh;
using MGroup.XFEM.Output.Vtk;

namespace MGroup.XFEM.FEM.Output
{
    public class ContinuousOutputMesh : IOutputMesh
    {
        private readonly IReadOnlyList<ICell<Node>> originalCells;
        private readonly IReadOnlyList<Node> originalVertices;
        private readonly List<VtkCell> outCells;
        private readonly List<VtkPoint> outVertices;

        public ContinuousOutputMesh(IReadOnlyList<Node> originalVertices, IReadOnlyList<ICell<Node>> originalCells)
        {
            this.originalVertices = originalVertices;
            this.originalCells = originalCells;

            var original2OutVertices = new Dictionary<Node, VtkPoint>();

            this.outVertices = new List<VtkPoint>();
            foreach (Node vertex in originalVertices)
            {
                var outVertex = new VtkPoint(vertex.ID, vertex.Coordinates);
                outVertices.Add(outVertex);
                original2OutVertices[vertex] = outVertex;
            }

            this.outCells = new List<VtkCell>();
            foreach (ICell<Node> cell in originalCells)
            {
                List<VtkPoint> vertices = cell.Nodes.Select(v => original2OutVertices[v]).ToList();
                outCells.Add(new VtkCell(cell.CellType, vertices));
            }
        }

        public int NumOutCells => outCells.Count;

        public int NumOutVertices => outVertices.Count;

        public IEnumerable<ICell<Node>> OriginalCells => originalCells;
        public IEnumerable<ICell<Node>> OriginalCellsList => originalCells;

        /// <summary>
        /// Same order as the corresponding one in <see cref="OutVertices"/>.
        /// </summary>
        public IEnumerable<Node> OriginalVertices => originalVertices;
        public IReadOnlyList<Node> OriginalVerticesList => originalVertices;

        public IEnumerable<VtkCell> OutCells => outCells;
        public List<VtkCell> CellsList => outCells;

        /// <summary>
        /// Same order as the corresponding one in <see cref="OriginalVertices"/>.
        /// </summary>
        public IEnumerable<VtkPoint> OutVertices => outVertices;
        public List<VtkPoint> VerticesList => outVertices;
    }
}
