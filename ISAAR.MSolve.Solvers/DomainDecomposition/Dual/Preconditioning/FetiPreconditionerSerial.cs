using System;
using System.Collections.Generic;
using System.Text;
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
    public class FetiPreconditionerSerial : IFetiPreconditioner
    {
        private readonly IFetiMatrixManager matrixManager;
        private readonly Dictionary<ISubdomain, IMappingMatrix> matricesBpb;
        private readonly IModel model;
        private readonly IFetiPreconditioningOperations preconditioning;

        private FetiPreconditionerSerial(IFetiPreconditioningOperations preconditioning, IModel model, 
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
            IFetiMatrixManager matrixManager, IStiffnessDistribution stiffnessDistribution)
        {
            this.preconditioning = preconditioning;
            this.model = model;
            this.matrixManager = matrixManager;

            this.matricesBpb = new Dictionary<ISubdomain, IMappingMatrix>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                matricesBpb[subdomain] = PreconditioningUtilities.CalcBoundaryPreconditioningBooleanMatrix(
                    subdomain, dofSeparator, lagrangeEnumerator, stiffnessDistribution); //TODO: When can these ones be reused?
                preconditioning.PrepareSubdomainSubmatrices(matrixManager.GetSubdomainMatrixManager(subdomain));
            }
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            foreach (ISubdomain sub in matricesBpb.Keys)
            {
                IFetiSubdomainMatrixManager subdomainMatrices = matrixManager.GetSubdomainMatrixManager(sub);
                Vector subdomainContribution = 
                    preconditioning.PreconditionSubdomainVector(rhs, subdomainMatrices, matricesBpb[sub]);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs) // bookmark1
        {
            foreach (ISubdomain sub in matricesBpb.Keys)
            {
                IFetiSubdomainMatrixManager subdomainMatrices = matrixManager.GetSubdomainMatrixManager(sub);
                Matrix subdomainContribution =
                    preconditioning.PreconditionSubdomainMatrix(rhs, subdomainMatrices, matricesBpb[sub]);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public class Factory : IFetiPreconditionerFactory
        {
            public IFetiPreconditioner CreatePreconditioner(IFetiPreconditioningOperations preconditioning,
                IModel model, IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                IFetiMatrixManager matrixManager, IStiffnessDistribution stiffnessDistribution)
            {
                return new FetiPreconditionerSerial(preconditioning, model, dofSeparator, lagrangeEnumerator, matrixManager, 
                    stiffnessDistribution);
            }
        }
    }
}
