using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public interface IAugmentationConstraintsFactory
    {
        IAugmentationConstraints CreateAugmentationConstraints(IModel model, IMidsideNodesSelection midsideNodesSelection,
            IFetiDPDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangesEnumerator);
    }
}
