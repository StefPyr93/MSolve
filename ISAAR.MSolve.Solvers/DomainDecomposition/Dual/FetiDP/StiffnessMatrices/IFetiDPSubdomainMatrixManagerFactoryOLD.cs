using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public interface IFetiDPSubdomainMatrixManagerFactoryOLD
    {
        //TODO: This is in a different namespace
        IFetiDPCoarseProblemSolverOLD CreateCoarseProblemSolver(IModel model); 

        IFetiDPSubdomainMatrixManagerOLD CreateMatricesManager(ISubdomain subdomain);
    }
}
