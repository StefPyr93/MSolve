using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This should be an enum class. There are only 2 possible cases.
//TODO: This should work for both FETI-1 and FETI-DP
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution 
{
    public interface IStiffnessDistribution : INodalLoadDistributor
    {
        double[] CalcBoundaryDofCoefficients(ISubdomain subdomain);

        IMappingMatrix CalcBoundaryPreconditioningSignedBooleanMatrix(ILagrangeMultipliersEnumerator lagrangeEnumerator, 
            ISubdomain subdomain, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix); //TODO: LagrangeEnumerator is only useful for heterogeneous. It should be injected in that contructor.

        void Update();
    }
}
