using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public interface IFetiSolver : ISolverMpi
    {
        Vector previousLambda { get; set; }
        bool usePreviousLambda { get; set; }

        IFetiDPInterfaceProblemSolver InterfaceProblemSolver { get; }
    }
}
