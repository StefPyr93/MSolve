using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public class LumpedPreconditionerOLD : IFetiPreconditioner
    {
        private readonly Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagers;
        private readonly Dictionary<int, IMappingMatrix> preconditioningBoundarySignedBooleanMatrices;
        private readonly int[] subdomainIDs;

        private LumpedPreconditionerOLD(int[] subdomainIDs, Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagers,
            Dictionary<int, IMappingMatrix> preconditioningBoundarySignedBooleanMatrices)
        {
            this.subdomainIDs = subdomainIDs;
            this.matrixManagers = matrixManagers;
            this.preconditioningBoundarySignedBooleanMatrices = preconditioningBoundarySignedBooleanMatrices;
        }

        //TODO: This can be moved to a base class. Only the S matrix is different for these preconditioners. 
        //      Other ones might be different though.
        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            foreach (int s in subdomainIDs)
            {
                IFetiSubdomainMatrixManagerOLD matrixManager = matrixManagers[s];
                IMappingMatrix Bpb = preconditioningBoundarySignedBooleanMatrices[s];

                // inv(F) * y = Bpb * Kbb * Bpb^T * y
                Vector temp = Bpb.Multiply(rhs, true);
                temp = matrixManager.MultiplyKbbTimes(temp);
                Vector subdomainContribution = Bpb.Multiply(temp);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            foreach (int s in subdomainIDs)
            {
                IFetiSubdomainMatrixManagerOLD matrixManager = matrixManagers[s];
                IMappingMatrix Bpb = preconditioningBoundarySignedBooleanMatrices[s];

                // inv(F) * y: Bpb * Kbb * Bpb^T * Y
                Matrix temp = Bpb.MultiplyRight(rhs, true);
                temp = matrixManager.MultiplyKbbTimes(temp);
                Matrix subdomainContribution = Bpb.MultiplyRight(temp);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public class Factory : FetiPreconditionerFactoryBaseOLD
        {
            public override bool ReorderInternalDofsForFactorization => false;

            public override IFetiPreconditioner CreatePreconditioner(IModel model, 
                IStiffnessDistributionOLD stiffnessDistribution, IDofSeparator dofSeparator, 
                ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator, Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagers)
            {
                int[] subdomainIDs = matrixManagers.Keys.ToArray();
                Dictionary<int, IMappingMatrix> boundaryBooleans = CalcBoundaryPreconditioningBooleanMatrices(model,
                    stiffnessDistribution, dofSeparator, lagrangeEnumerator);

                foreach (int s in subdomainIDs)
                {
                    ISubdomain subdomain = model.GetSubdomain(s);
                    if (!subdomain.StiffnessModified) continue;
                    Debug.WriteLine($"{typeof(DiagonalDirichletPreconditionerOLD).Name}.{this.GetType().Name}:"
                        + $" Extracting boundary/internal submatrices of subdomain {s} for preconditioning");
                    int[] boundaryDofs = dofSeparator.GetBoundaryDofIndices(subdomain);
                    matrixManagers[s].ExtractKbb(boundaryDofs);
                }

                return new LumpedPreconditionerOLD(subdomainIDs, matrixManagers, boundaryBooleans);
            }
        }
    }
}
