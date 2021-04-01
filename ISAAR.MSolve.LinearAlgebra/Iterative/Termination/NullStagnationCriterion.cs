using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    public class NullStagnationCriterion : IStagnationCriterion
    {
        public bool HasStagnated() => false;

        public void StoreInitialError(double initialError) { }

        public void StoreNewError(double currentError) { }

    }
}
