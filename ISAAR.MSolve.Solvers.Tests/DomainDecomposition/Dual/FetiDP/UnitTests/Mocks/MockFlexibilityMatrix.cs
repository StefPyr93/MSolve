using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks
{
    public class MockFlexibilityMatrix : IFetiDPFlexibilityMatrix
    {
        public int NumGlobalLagrangeMultipliers => Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers.Length;

        public Vector MultiplyFIrc(Vector vIn)
        {
            if (vIn != null) return Example4x4QuadsHomogeneous.MatrixFIrc * vIn;
            else return null;
        }

        public Vector MultiplyFIrcTransposed(Vector vIn)
        {
            if (vIn != null) return Example4x4QuadsHomogeneous.MatrixFIrc.Multiply(vIn, true);
            else return null;
        }

        public void MultiplyFIrr(Vector vIn, Vector vOut)
        {
            if ((vIn != null) && (vOut != null)) vOut.CopyFrom(Example4x4QuadsHomogeneous.MatrixFIrr * vIn);
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vIn)
        {
            if (vIn != null)
            {
                Vector FIrrTimesVector = Example4x4QuadsHomogeneous.MatrixFIrr * vIn;
                Vector FIrcTransposedTimesVector = Example4x4QuadsHomogeneous.MatrixFIrc.Multiply(vIn, true);
                return (FIrrTimesVector, FIrcTransposedTimesVector);
            }
            else return (null, null);
        }
    }
}
