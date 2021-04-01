using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;

//TODO: In XFEM where there might be dynamically created dofs (e.g. due to branching cracks), there should be a way to broadcast 
//      the dof serializer.
namespace ISAAR.MSolve.Discretization.Transfer
{
    public interface IDofSerializer
    {
        IDofType Deserialize(int dofID);
        int Serialize(IDofType dofType);
    }
}
