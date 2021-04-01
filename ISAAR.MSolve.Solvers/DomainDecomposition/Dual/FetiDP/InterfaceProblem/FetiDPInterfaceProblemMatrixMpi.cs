using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;

//TODO: perhaps I PCG should work with delegates instead of ILinearTransformation and FetiDPInterfaceProblemMatrixMpi/Serial 
//      should be a method in FetiDPInterfaceProblemSolverMpi/Serial.
//TODO: Remove duplication between this and the serial implementation. It may not be possible here though, since global operations
//      that are meant to be run only on master and global operations that are meant to be run everywhere are intermingled.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// This class concerns global operations and uses global data. In an MPI environment, instances of it should not be used in 
    /// processes other than master. 
    /// </summary>
    public class FetiDPInterfaceProblemMatrixMpi : ILinearTransformation
    {
        private readonly IFetiDPFlexibilityMatrix flexibility;
        private readonly IFetiDPMatrixManager matrixManager;
        private readonly ProcessDistribution procs;

        public FetiDPInterfaceProblemMatrixMpi(ProcessDistribution processDistribution, IFetiDPMatrixManager matrixManager, 
            IFetiDPFlexibilityMatrix flexibility)
        {
            this.procs = processDistribution;
            this.matrixManager = matrixManager;
            this.flexibility = flexibility;
        }

        public int NumColumns => flexibility.NumGlobalLagrangeMultipliers;
        public int NumRows => flexibility.NumGlobalLagrangeMultipliers;

        public void Multiply(IVectorView lhsVector, IVector rhsVector)
        {
            //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
            //      IVectorView and IVector
            var x = (Vector)lhsVector;
            var y = (Vector)rhsVector;

            // y = (FIrr + FIrc * inv(KccStar) * FIrc^T) * x
            (Vector FIrrX, Vector transpFIrcTranspX) = flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(x);
            Vector temp = null;
            if (procs.IsMasterProcess)
            {
                y.CopyFrom(FIrrX);
                temp = matrixManager.MultiplyInverseCoarseProblemMatrix(transpFIrcTranspX);
            }
            temp = flexibility.MultiplyFIrc(temp);
            if (procs.IsMasterProcess) y.AddIntoThis(temp);
        }
    }
}
