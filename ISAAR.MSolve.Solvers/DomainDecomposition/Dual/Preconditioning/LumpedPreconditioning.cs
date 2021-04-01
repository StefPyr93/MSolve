using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public class LumpedPreconditioning : IFetiPreconditioningOperations
    {
        public bool ReorderInternalDofsForFactorization => false;

        public void PrepareSubdomainSubmatrices(IFetiSubdomainMatrixManager matrixManager)
        {
            ISubdomain subdomain = matrixManager.LinearSystem.Subdomain;
            if (subdomain.StiffnessModified)
            {
                Debug.WriteLine($"{this.GetType().Name}:"
                    + $" Extracting boundary/internal submatrices of subdomain {subdomain.ID} for preconditioning");
                matrixManager.ExtractKbb();
            }
        }

        public Matrix PreconditionSubdomainMatrix(Matrix rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb)
        {
            // inv(F) * Y: Bpb * Kbb * Bpb^T * Y
            Matrix temp = Bpb.MultiplyRight(rhs, true);
            temp = matrixManager.MultiplyKbbTimes(temp);
            Matrix subdomainContribution = Bpb.MultiplyRight(temp);
            return subdomainContribution;
        }

        public Vector PreconditionSubdomainVector(Vector rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb)
        {
            // inv(F) * y = Bpb * Kbb * Bpb^T * y
            Vector temp = Bpb.Multiply(rhs, true);
            temp = matrixManager.MultiplyKbbTimes(temp);
            Vector subdomainContribution = Bpb.Multiply(temp);
            return subdomainContribution;
        }
    }
}
