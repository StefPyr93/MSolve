using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INodalLoad
    {
        INode Node { get; }
        IDofType DOF { get; } //TODO: Rename to Dof
        double Amount { get; }
    }
}
