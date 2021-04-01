using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public class FetiDPFlexibilityMatrixSerial : IFetiDPFlexibilityMatrix
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix> subdomainFlexibilities;

        public FetiDPFlexibilityMatrixSerial(IModel model, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPMatrixManager matrixManager)
        {
            this.dofSeparator = dofSeparator;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.NumGlobalLagrangeMultipliers = lagrangesEnumerator.NumLagrangeMultipliers;

            this.subdomainFlexibilities = new Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                this.subdomainFlexibilities[sub] = new FetiDPSubdomainFlexibilityMatrix(sub, dofSeparator, lagrangesEnumerator,
                    matrixManager);
            }
        }

        public int NumGlobalLagrangeMultipliers { get; }

        public Vector MultiplyFIrc(Vector vIn)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrc(vIn, dofSeparator);
            var vOut = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrc(vIn);
                vOut.AddIntoThis(subdomainRhs); //TODO: The operation y = A * x + y is important here to avoid extra allocations and additions
            }
            return vOut;
        }

        public Vector MultiplyFIrcTransposed(Vector vIn)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            var vOut = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrcTransposed(vIn);
                vOut.AddIntoThis(subdomainRhs); //TODO: The operation y = A * x + y is important here to avoid extra allocations and additions
            }
            return vOut;
        }

        public void MultiplyFIrr(Vector vIn, Vector vOut)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrr(vIn, vOut, lagrangesEnumerator);
            foreach (ISubdomain sub in subdomainFlexibilities.Keys)
            {
                Vector subdomainRhs = subdomainFlexibilities[sub].MultiplyFIrr(vIn);
                vOut.AddIntoThis(subdomainRhs); //TODO: The operation y = A * x + y is important here to avoid extra allocations and additions
            }
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vIn)
        {
            FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            var FIrrTimesVector = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            var FIrcTransposedTimesVector = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs);
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