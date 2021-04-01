using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

//TODO: Also add ReorderInternalDofsForMultiplication and ReorderBoundaryDofsForMultiplication
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public interface IFetiPreconditionerFactoryOLD
    {
        bool ReorderInternalDofsForFactorization { get; }

        IFetiPreconditioner CreatePreconditioner(IModel model, IStiffnessDistributionOLD stiffnessDistribution,
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator, 
            Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagers);
    }
}
