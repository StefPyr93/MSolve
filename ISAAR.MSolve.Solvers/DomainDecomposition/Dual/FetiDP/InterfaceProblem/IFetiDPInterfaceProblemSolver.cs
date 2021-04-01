using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;

//TODO: It would be better if globalForcesNorm was read by implementations of this interface from another object. Right now,
//      in an MPI environment only master process has a value to pass. All others can pass some dummy value that will be ignored.
//      The same holds for the returned values.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPInterfaceProblemSolver
    {
        Vector PreviousLambda { get; set; }
        bool UsePreviousLambda { get; set; }
        ReorthogonalizedPcg Pcg { get; }

        bool UseStagnationCriterion { get; set; }

        Vector SolveInterfaceProblem(IFetiDPMatrixManager matrixManager,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, double globalForcesNorm, ISolverLogger logger);
    }
}
