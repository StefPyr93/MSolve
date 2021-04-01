using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;

//TODO: This could be split into an interface with the same name and an IFtiDPCoarseProblemSolver.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPInterfaceProblemSolverOLD
    {
        (Vector lagrangeMultipliers, Vector cornerDisplacements) SolveInterfaceProblem(FetiDPFlexibilityMatrixOLD flexibility, 
            IFetiPreconditioner preconditioner, IFetiDPCoarseProblemSolverOLD coarseProblemSolver, Vector globalFcStar, Vector dr,
            double globalForcesNorm, SolverLoggerOLD logger);
    }
}
