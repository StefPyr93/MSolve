using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices
{
    public class FetiDP3dMatrixManagerFactoryDense : IFetiDP3dMatrixManagerFactory
    {
        public IFetiDPGlobalMatrixManager CreateGlobalMatrixManager(IModel model, IFetiDPDofSeparator dofSeparator,
            IAugmentationConstraints augmentationConstraints) 
            => new FetiDP3dGlobalMatrixManagerDense(model, dofSeparator, augmentationConstraints);

        public IFetiDPSubdomainMatrixManager CreateSubdomainMatrixManager(ISubdomain subdomain, 
            IFetiDPDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangesEnumerator, 
            IAugmentationConstraints augmentationConstraints)
            => new FetiDP3dSubdomainMatrixManagerDense(subdomain, dofSeparator, lagrangesEnumerator, augmentationConstraints);
    }
}
