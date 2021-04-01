using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;

//TODO: Combine this with convergence criterion
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    public interface IStagnationCriterion
    {
        bool HasStagnated();

        void StoreInitialError(double initialError);

        void StoreNewError(double currentError);

    }
}
