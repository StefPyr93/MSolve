using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: decide on consistent collections, indexes/keys and naming
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    public interface ICrackDescriptionMpi : ICrackDescription
    {
        void ScatterCrackData(IXModelMpi model);
    }
}
