using System;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Transfer
{
    [Serializable]
    public class XElementDto
    {
        public CellType cellType;
        public int id;
        public int[] nodeIds;

        public XElementDto(IXFiniteElement element)
        {
            this.id = element.ID;
            this.cellType = ((IElementType)element).CellType;
            this.nodeIds = new int[element.Nodes.Count];
            for (int n = 0; n < nodeIds.Length; ++n) nodeIds[n] = element.Nodes[n].ID;
        }

        public IXFiniteElement Deserialize(IXFiniteElementFactory elementFactory, Func<int, XNode> getNodeWithID)
        {
            var nodes = new XNode[nodeIds.Length];
            for (int n = 0; n < nodeIds.Length; ++n) nodes[n] = getNodeWithID(nodeIds[n]);
            return elementFactory.CreateElement(id, cellType, nodes);
        }
    }
}
