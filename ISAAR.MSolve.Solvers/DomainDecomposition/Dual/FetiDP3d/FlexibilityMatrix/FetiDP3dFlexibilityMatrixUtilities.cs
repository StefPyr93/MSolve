using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This should not exist. Its code should be defined in FetiDPFlexibilityMatrixBase. It is not like other CPW Part classes
//      which actually stored state.
//TODO: Useless checks probably. Should be removed
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix
{
    public static class FetiDP3dFlexibilityMatrixUtilities
    {
        [Conditional("DEBUG")]
        public static void CheckMultiplicationGlobalFIrc3d(Vector vIn, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator, IAugmentationConstraints augmentationConstraints)
        {
            Preconditions.CheckMultiplicationDimensions(
                dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints, vIn.Length);
        }

        [Conditional("DEBUG")]
        public static void CheckMultiplicationGlobalFIrcTransposed(Vector vIn, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangeEnumerator)
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
