using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

//TODO: Perhaps I do not need factories with the new design. I could call directly new FetiPreconditionerMpi(...) in the solver. 
//      IFetiPreconditioningOperations are also very easy for the user. The factories would be useful if there was only 1 solver 
//      class that transparently used MPI/serial implementations of its components, but that is too optimistic.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public interface IFetiPreconditionerFactory
    {
        IFetiPreconditioner CreatePreconditioner(IFetiPreconditioningOperations preconditioning,
            IModel model, IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
            IFetiMatrixManager matrixManager, IStiffnessDistribution stiffnessDistribution);
    }
}
