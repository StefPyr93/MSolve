using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPMatrixManagerFactoryDense : IFetiDPMatrixManagerFactory
    {
        public IFetiDPGlobalMatrixManager CreateGlobalMatrixManager(IModel model, IFetiDPDofSeparator dofSeparator) 
            => new FetiDPGlobalMatrixManagerDense(model, dofSeparator);

        public IFetiDPSubdomainMatrixManager CreateSubdomainMatrixManager(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator)
            => new FetiDPSubdomainMatrixManagerDense(subdomain, dofSeparator);
    }
}
