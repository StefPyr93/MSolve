using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    internal static class PreconditioningUtilities
    {
        //TODO: Is it worth having a whole class for these 3 lines?
        internal static IMappingMatrix CalcBoundaryPreconditioningBooleanMatrix(ISubdomain subdomain, IDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator, IStiffnessDistribution stiffnessDistribution)
        {
            SignedBooleanMatrixColMajor B = lagrangeEnumerator.GetBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Bb = B.GetColumns(dofSeparator.GetBoundaryDofIndices(subdomain), false);
            return stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrix(lagrangeEnumerator, subdomain, Bb);
        }
    }
}
