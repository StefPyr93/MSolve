using ISAAR.MSolve.Analyzers.Loading;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzerProvider
    {
        IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        void ClearMatrices();
        void Reset();
    }
}
