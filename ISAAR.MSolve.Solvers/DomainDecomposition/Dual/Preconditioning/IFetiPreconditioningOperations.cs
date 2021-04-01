using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public interface IFetiPreconditioningOperations
    {
        bool ReorderInternalDofsForFactorization { get; }

        Matrix PreconditionSubdomainMatrix(Matrix rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb);
        Vector PreconditionSubdomainVector(Vector rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb);

        void PrepareSubdomainSubmatrices(IFetiSubdomainMatrixManager matrixManager);
    }
}
