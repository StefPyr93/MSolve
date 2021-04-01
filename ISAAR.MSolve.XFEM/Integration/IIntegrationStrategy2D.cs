using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Materials;

//TODO: Not sure if generics are needed. I could just have different interfaces for different element types. 
//      It would make sense if some of the actual integration strategy implementations would work for many element types.
//TODO: Actually the generics create complications when I need to transfer the integration rules in an MPI environment.
namespace ISAAR.MSolve.XFEM.Integration
{
    /// <summary>
    /// Algorithms for complex integration rules for specific finite element types. These need the data from each 
    /// finite element to generate integration points for use only by that finite element. 
    /// They typically make use of the standard quadrature rules.
    /// </summary>
    /// <typeparam name="TElement"></typeparam>
    public interface IIntegrationStrategy2D<TElement>
    {
        IReadOnlyList<GaussPoint> GenerateIntegrationPoints(TElement element);
    }
}
