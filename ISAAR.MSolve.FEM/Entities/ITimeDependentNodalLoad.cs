using ISAAR.MSolve.Discretization.FreedomDegrees;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: The design of time dependent nodal loads should be as follows:
//      - Name them ITransientNodalLoad
//      - They implement the same interface as steady state nodal loads. This is very important.
//          - The same would apply for non nodal loads.
//      - They hold a reference to an ITimeFunction that has a method Evaluate().
//      - They return Amount = baseAmount * timeFunc.Evaluate()
//      - ITimeFunction can be shared among many loads, so that it only needs to be evaluated once. This means it stores state.
//      - ITimeFunction is an observer of ITransientModel. 
//      - ITransientModel is updated by the analyzer and it notifies its observers to update their state.
namespace ISAAR.MSolve.FEM.Entities
{
    public interface ITimeDependentNodalLoad
    {
        Node Node { get; set; }
        IDofType DOF { get; set; }

        double GetLoadAmount(int timeStep);
    }
}
