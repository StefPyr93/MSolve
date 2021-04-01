using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Transfer.Elements
{
    public interface IElementDto
    {
        Element Deserialize(IReadOnlyDictionary<int, Node> allNodes, Dictionary<int, IFiniteElementMaterial> allMaterials);
    }
}
