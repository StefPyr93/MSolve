using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;

//TODO: This class should not exist. The FetiDPFlexibilityMatrixSerial used in regular FETI-DP should be used instead. However 
//      rename it to match the more general 3D names
//TODO: The serial/MPI coordinators of regular FETI-DP should be used instead of this. Only the subdomain operations and 
//      the dimensions are different.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix
{
    public class FetiDP3dFlexibilityMatrixSerial : IFetiDPFlexibilityMatrix
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix> subdomainFlexibilities;

        public FetiDP3dFlexibilityMatrixSerial(IModel model, IFetiDPDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IAugmentationConstraints augmentationConstraints, 
            IFetiDPMatrixManager matrixManager) 
        {
            this.dofSeparator = dofSeparator;
            this.augmentationConstraints = augmentationConstraints;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.NumGlobalLagrangeMultipliers = lagrangesEnumerator.NumLagrangeMultipliers;

            this.subdomainFlexibilities = new Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                this.subdomainFlexibilities[sub] = new FetiDP3dSubdomainFlexibilityMatrix(sub, dofSeparator, lagrangesEnumerator,
                    augmentationConstraints, matrixManager);
            }
        }

        public int NumGlobalLagrangeMultipliers { get; }

        public Vector MultiplyFIrc(Vector vIn)
        {
            FetiDP3dFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrc3d(vIn, dofSeparator, lagrangesEnumerator, 
                augmentationConstraints);
            var vOut = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrc(vIn);
                vOut.AddIntoThis(subdomainRhs);
            }
            return vOut;
        }

        public Vector MultiplyFIrcTransposed(Vector vIn)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            var vOut = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrcTransposed(vIn);
                vOut.AddIntoThis(subdomainRhs);
            }
            return vOut;
        }

        public void MultiplyFIrr(Vector vIn, Vector vOut)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrr(vIn, vOut, lagrangesEnumerator);
            vOut.Clear();
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrr(vIn);
                vOut.AddIntoThis(subdomainRhs);
            }
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vIn)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            var FIrrTimesVector = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            var FIrcTransposedTimesVector = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                (Vector subdomainFIrrTimesVector, Vector subdomainFIrcTransposedTimesVector) =
                    subdomainFlexibilities[sub].MultiplyFIrrAndFIrcTransposedTimesVector(vIn);
                FIrrTimesVector.AddIntoThis(subdomainFIrrTimesVector); //TODO: The operation y = A * x + y is important here to avoid extra allocations and additions
                FIrcTransposedTimesVector.AddIntoThis(subdomainFIrcTransposedTimesVector);
            }
            return (FIrrTimesVector, FIrcTransposedTimesVector);
        }
    }
}
