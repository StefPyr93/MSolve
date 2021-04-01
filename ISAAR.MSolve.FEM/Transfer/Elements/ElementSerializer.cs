using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: going through all the if clauses for each element is slow. Use a dictionary instead.
namespace ISAAR.MSolve.FEM.Transfer.Elements
{
    public class ElementSerializer
    {
        public IElementDto Serialize(Element element)
        {
            if (element.ElementType is ContinuumElement2D continuum)
            {
                if (continuum.Interpolation is InterpolationQuad4) return new Quad4Dto(element);
            }
            throw new ArgumentException("Unknown element type.");
        }
    }
}
