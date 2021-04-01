using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPMatrixManagerFactorySuitesparse : IFetiDPMatrixManagerFactory
    {
        private readonly IReorderingAlgorithm reordering;

        public FetiDPMatrixManagerFactorySuitesparse()
        {
            this.reordering = new OrderingAmdSuiteSparse();
        }

        public FetiDPMatrixManagerFactorySuitesparse(IReorderingAlgorithm reordering)
        {
            this.reordering = reordering;
        }

        public IFetiDPGlobalMatrixManager CreateGlobalMatrixManager(IModel model, IFetiDPDofSeparator dofSeparator) 
            => new FetiDPGlobalMatrixManagerSuiteSparse(model, dofSeparator, reordering);

        public IFetiDPSubdomainMatrixManager CreateSubdomainMatrixManager(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator)
            => new FetiDPSubdomainMatrixManagerSuiteSparse(subdomain, dofSeparator, reordering);
    }
}
