using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public interface ILagrangeMultipliersEnumerator
    {
        IReadOnlyList<LagrangeMultiplier> LagrangeMultipliers { get; }

        int NumLagrangeMultipliers { get; }

        ///// <summary>
        ///// Associates each lagrange multiplier with the instances of the boundary dof, for which continuity is enforced. 
        ///// </summary>
        //IEnumerable<LagrangeMultiplier> EnumerateOrderedLagrangeMultipliers();

        SignedBooleanMatrixColMajor GetBooleanMatrix(ISubdomain subdomain);
    }
}
