using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rename to FetiDPCoarseProblemSolver
//TODO: Not sure about the interface. These methods should be wrapped by IFetiDPMatrixManager
//TODO: Calculating the coarse problem rhs is also subdject to mpi/serial implementations, but it does not depend on any matrix format.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public interface IFetiDPGlobalMatrixManager 
    {
        Vector CoarseProblemRhs { get; }

        void CalcCoarseProblemRhs(Dictionary<ISubdomain, Vector> condensedRhsVectors);

        //TODO: Does not make sense to only provide cornerNodeSelection, without midsideNodeSelection. However the latter is 
        //      injected through the constructor
        void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection,
            Dictionary<ISubdomain, IMatrixView> subdomainCoarseMatrices);

        void ClearCoarseProblemRhs();
        void ClearInverseCoarseProblemMatrix();

        Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector);

        DofPermutation ReorderGlobalCornerDofs();
    }
}
