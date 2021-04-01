using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This class should be safe to call methods from, regardless which process it is.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public interface IFetiDPFlexibilityMatrix
    {
        int NumGlobalLagrangeMultipliers { get; }

        Vector MultiplyFIrc(Vector vIn);
        Vector MultiplyFIrcTransposed(Vector vIn);
        void MultiplyFIrr(Vector vIn, Vector vOut); //TODO: Unused, contrary to FIrc^T

        (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vIn);

    }
}
