using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This should be an enum class. There are only 2 possible cases.
//TODO: This should work for both FETI-1 and FETI-DP
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution 
{
    public interface IStiffnessDistributionOLD : INodalLoadDistributor
    {
        double[] CalcBoundaryDofCoefficients(ISubdomain subdomain);

        Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator, 
            Dictionary<int, SignedBooleanMatrixColMajor> boundarySignedBooleanMatrices);

        void Update(Dictionary<int, IMatrixView> stiffnessesFreeFree);

        //Dictionary<int, SparseVector> DistributeNodalLoadsOLD(Dictionary<int, ISubdomain> subdomains,
        //    Table<INode, IDofType, double> globalNodalLoads);
    }
}
