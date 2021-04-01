using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: going through all the if clauses for each material is slow. Use a dictionary instead.
namespace ISAAR.MSolve.Materials.Transfer
{
    public class MaterialSerializer
    {
        public IMaterialDto Serialize(IFiniteElementMaterial material)
        {
            if (material is ElasticMaterial2D elastic2D) return new Elastic2DMaterialDto(elastic2D);
            throw new ArgumentException("Unknown material type.");
        }
    }
}
