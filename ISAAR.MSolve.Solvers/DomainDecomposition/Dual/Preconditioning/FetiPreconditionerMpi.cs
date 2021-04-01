using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.LinearAlgebra.Distributed.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public class FetiPreconditionerMpi : IFetiPreconditioner
    {
        private readonly Dictionary<int, IMappingMatrix> matricesBpb;
        private readonly Dictionary<int, IFetiSubdomainMatrixManager> matrixManagers;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiPreconditioningOperations preconditioning;
        private readonly ProcessDistribution procs;

        private FetiPreconditionerMpi(ProcessDistribution processDistribution, IFetiPreconditioningOperations preconditioning,
            IModel model, IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangesEnumerator,
            IFetiMatrixManager matrixManager, IStiffnessDistribution stiffnessDistribution)
        {
            this.procs = processDistribution;
            this.preconditioning = preconditioning;
            this.lagrangesEnumerator = lagrangesEnumerator;

            this.matricesBpb = new Dictionary<int, IMappingMatrix>();
            this.matrixManagers = new Dictionary<int, IFetiSubdomainMatrixManager>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                this.matricesBpb[s] = PreconditioningUtilities.CalcBoundaryPreconditioningBooleanMatrix(
                    subdomain, dofSeparator, lagrangesEnumerator, stiffnessDistribution); //TODO: When can these ones be reused?
                this.matrixManagers[s] = matrixManager.GetSubdomainMatrixManager(subdomain);
                preconditioning.PrepareSubdomainSubmatrices(matrixManager.GetSubdomainMatrixManager(subdomain));
            }
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref rhs, lagrangesEnumerator.NumLagrangeMultipliers);
            var subdomainContributions = new List<Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                subdomainContributions.Add(preconditioning.PreconditionSubdomainVector(rhs, matrixManagers[s], matricesBpb[s]));
            }
            transferrer.SumVectors(subdomainContributions, lhs);
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs)
        {
            var transferrer = new MatrixTransferrer(procs);
            transferrer.BroadcastMatrix(ref rhs); //TODO: Perhaps the dimensions are already known and we can avoid broadcasting them too.
            var subdomainContributions = new List<Matrix>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                subdomainContributions.Add(preconditioning.PreconditionSubdomainMatrix(rhs, matrixManagers[s], matricesBpb[s]));
            }
            transferrer.SumMatrices(subdomainContributions, lhs);
        }

        public class Factory : IFetiPreconditionerFactory
        {
            private readonly ProcessDistribution procs;

            public Factory(ProcessDistribution processDistribution)
            {
                this.procs = processDistribution;
            }

            public IFetiPreconditioner CreatePreconditioner(IFetiPreconditioningOperations preconditioning,
                IModel model, IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                IFetiMatrixManager matrixManager, IStiffnessDistribution stiffnessDistribution)
            {
                return new FetiPreconditionerMpi(procs, preconditioning, model, dofSeparator, lagrangeEnumerator, matrixManager, 
                    stiffnessDistribution);
            }
        }
    }
}
