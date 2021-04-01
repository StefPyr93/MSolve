using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    ////TODO: Using these 2 and the methods that take them as arguments, I could write classes that have nothing to do with the DDM
    ////      itself. They only know how to assemble matrices and vectors in serial/MPI environment. I can take it even further
    ////      and have the SubdomainMatrixAssembler use the serial implementation.
    //public delegate (IMatrixView subdomainMatrix, UnsignedBooleanMatrix globalToSubdomainDofsMap) CalcSubdomainMatrix(
    //    ISubdomain subdomain);
    //public delegate (Vector subdomainVector, UnsignedBooleanMatrix globalToSubdomainDofsMap) CalcSubdomainVector(
    //    ISubdomain subdomain);

    public interface IFetiDPMatrixManager : IFetiMatrixManager, IFetiDPSeparatedDofReordering
    {

        //TODO: Not sure about this. I should probably wrap the methods in this interfaces, in order to avoid having them 
        //      called in a loop. E.g MultiplyMatrixTimes... of preconditioner. In this case the in and out vectors should also 
        //      be contained in classes with MPI/Serial implementation or at least in LinearSystem-like classes, so that each 
        //      process read and writes to its corresponding one.
        IFetiDPSubdomainMatrixManager GetFetiDPSubdomainMatrixManager(ISubdomain subdomain);

        Vector CoarseProblemRhs { get; }

        void CalcCoarseProblemRhs();
        void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection);

        void ClearCoarseProblemRhs();
        void ClearInverseCoarseProblemMatrix();

        Vector MultiplyInverseCoarseProblemMatrix(Vector vector);
    }
}
