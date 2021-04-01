using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    //TODO: Confusing name. The child analyzer of this is a nonlinear analyzer.
    public interface INonLinearParentAnalyzer : IParentAnalyzer
    {
        //TODO: Also confusing name: OtherRhsComponents. Other than what?
        IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution);
    }
}
