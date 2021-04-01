using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: perhaps IFeti1SubdomainMatrixManager and IFetiDPSubdomainMatrixManager should not inherit from this one.
//TODO: perhaps I should providing access to the assembler instead of wrapping its methods.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public interface IFetiSubdomainMatrixManager
    {
        ISingleSubdomainLinearSystemMpi LinearSystem { get; }

        void BuildFreeDofsMatrix(ISubdomainFreeDofOrdering dofOrdering, IElementMatrixProvider matrixProvider); //TODO: This could be called by the serial/MPI IFetiMatrixManager, instead of the solver.

        (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildFreeConstrainedMatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider);

        void ClearMatrices();

        void ExtractBoundaryInternalSubmatricesAndInvertKii(bool diagonalKii);

        void ExtractKbb();
        void HandleDofOrderingWillBeModified();

        Vector MultiplyInverseKiiTimes(Vector vector, bool diagonalOnly);
        Matrix MultiplyInverseKiiTimes(Matrix matrix, bool diagonalOnly);
        Vector MultiplyKbbTimes(Vector vector);
        Matrix MultiplyKbbTimes(Matrix matrix);
        Vector MultiplyKbiTimes(Vector vector);
        Matrix MultiplyKbiTimes(Matrix matrix);
        Vector MultiplyKibTimes(Vector vector);
        Matrix MultiplyKibTimes(Matrix matrix);
    }
}
