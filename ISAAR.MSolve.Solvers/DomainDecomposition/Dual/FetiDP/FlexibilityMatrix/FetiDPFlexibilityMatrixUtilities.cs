using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Useless checks probably. Should be removed. Either way it complicates things needlessly
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public static class FetiDPFlexibilityMatrixUtilities
    {
        [Conditional("DEBUG")]
        public static void CheckMultiplicationGlobalFIrc(Vector vIn, IFetiDPDofSeparator dofSeparator)
        {
            Preconditions.CheckMultiplicationDimensions(dofSeparator.NumGlobalCornerDofs, vIn.Length);
        }

        [Conditional("DEBUG")]
        public static void CheckMultiplicationGlobalFIrcTransposed(Vector vIn, ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            Preconditions.CheckMultiplicationDimensions(lagrangeEnumerator.NumLagrangeMultipliers, vIn.Length);
        }

        [Conditional("DEBUG")]
        public static void CheckMultiplicationGlobalFIrr(Vector vIn, Vector vOut, 
            ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            Preconditions.CheckMultiplicationDimensions(lagrangeEnumerator.NumLagrangeMultipliers, vIn.Length);
            Preconditions.CheckSystemSolutionDimensions(lagrangeEnumerator.NumLagrangeMultipliers, vOut.Length);
        }
    }
}
