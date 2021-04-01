using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.XFEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Elements
{
    public interface IXFiniteElementFactory
    {
        IXFiniteElement CreateElement(int id, CellType cellType, IReadOnlyList<XNode> nodes);
    }
}
