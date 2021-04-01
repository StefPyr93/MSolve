using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
//TODO: Reanalysis: if the global lagranges have not changed, Bpb can be reused.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public abstract class FetiPreconditionerFactoryBaseOLD : IFetiPreconditionerFactoryOLD
    {
        public abstract bool ReorderInternalDofsForFactorization { get; }

        public abstract IFetiPreconditioner CreatePreconditioner(IModel model, 
            IStiffnessDistributionOLD stiffnessDistribution, IDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator, Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagers);

        protected Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningBooleanMatrices(IModel model, 
            IStiffnessDistributionOLD stiffnessDistribution, IDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator)
        {
            var matricesBb = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[s];
                SignedBooleanMatrixColMajor Bb = B.GetColumns(dofSeparator.GetBoundaryDofIndices(subdomain), false);
                matricesBb[s] = Bb;
            }
            Dictionary<int, IMappingMatrix> matricesBpb = stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrices(
                lagrangeEnumerator, matricesBb);
            
            return matricesBpb;
        }
    }
}
