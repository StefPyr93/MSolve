using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Vectors;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: This class, interface should not exist. All these multiplications (e.g. FIrr * vector) should be handled by 
//      InterfaceProblemMatrix. (FIrr * FIrc * inv(KccStar) * FIrc^T) * vector should be a dedicated method, to remove the need 
//      for caching intermediate results.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public class FetiDPFlexibilityMatrixMpi : IFetiDPFlexibilityMatrix
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IModel model;
        private readonly ProcessDistribution procs;
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix> subdomainFlexibilities;

        public FetiDPFlexibilityMatrixMpi(ProcessDistribution procs, IModel model, IFetiDPDofSeparator dofSeparator,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPMatrixManager matrixManager)
        {
            this.procs = procs;
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.NumGlobalLagrangeMultipliers = lagrangesEnumerator.NumLagrangeMultipliers;

            this.subdomainFlexibilities = new Dictionary<ISubdomain, IFetiDPSubdomainFlexibilityMatrix>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                this.subdomainFlexibilities[subdomain] =
                    new FetiDPSubdomainFlexibilityMatrix(subdomain, dofSeparator, lagrangesEnumerator, matrixManager);
            }
        }

        public int NumGlobalLagrangeMultipliers { get; }

        /// <summary>
        /// The returned vector will be null in all processes other master.
        /// </summary>
        /// <param name="vIn">This vector must be present in all processes.</param>
        /// <returns></returns>
        public Vector MultiplyFIrc(Vector vIn)
        {
            if (procs.IsMasterProcess)
            {
                FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrc(vIn, dofSeparator);
            }
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref vIn, dofSeparator.NumGlobalCornerDofs);
            //IEnumerable<Vector> subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplyFIrc(vIn));
            Vector[] subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplyFIrc(vIn)).ToArray();
            return transferrer.SumVectors(subdomainRhs);
        }

        /// <summary>
        /// The returned vector will be null in all processes other master.
        /// </summary>
        /// <param name="vIn">This vector must be present in all processes.</param>
        /// <returns></returns>
        public Vector MultiplyFIrcTransposed(Vector vIn)
        {
            if (procs.IsMasterProcess)
            {
                FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            }
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref vIn, lagrangesEnumerator.NumLagrangeMultipliers);
            //IEnumerable<Vector> subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplySubdomainFIrcTransposed(vIn)); //TODO: This version calculates FIrc^T * vector twice (not FIrr thouh...)
            Vector[] subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplyFIrcTransposed(vIn)).ToArray();
            return transferrer.SumVectors(subdomainRhs);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vIn">This vector must be present in all processes.</param>
        /// <param name="vOut">It will be ignored in processes other than master.</param>
        public void MultiplyFIrr(Vector vIn, Vector vOut)
        {
            if (procs.IsMasterProcess)
            {
                FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrr(vIn, vOut, lagrangesEnumerator);
            }
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref vIn, lagrangesEnumerator.NumLagrangeMultipliers);
            //IEnumerable<Vector> subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplySubdomainFIrr(vIn)); //TODO: This version calculates FIrc^T * vector twice (not FIrr thouh...)
            Vector[] subdomainRhs = subdomainFlexibilities.Values.Select(F => F.MultiplyFIrr(vIn)).ToArray();
            transferrer.SumVectors(subdomainRhs, vOut);
        }

        public (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vIn)
        {
            if (procs.IsMasterProcess)
            {
                FetiDPFlexibilityMatrixUtilities.CheckMultiplicationGlobalFIrcTransposed(vIn, lagrangesEnumerator);
            }
            var transferrer = new VectorTransferrer(procs);
            transferrer.BroadcastVector(ref vIn, lagrangesEnumerator.NumLagrangeMultipliers);
            var FIrrTimesVectorPerSubdomain = new List<Vector>();
            var FIrcTransposedTimesVectorPerSubdomain = new List<Vector>();
            foreach (IFetiDPSubdomainFlexibilityMatrix flexibility in subdomainFlexibilities.Values)
            {
                (Vector FIrrTimesVectorOfSubdomain, Vector FIrcTransposedTimesVectorOfSubdomain) =
                    flexibility.MultiplyFIrrAndFIrcTransposedTimesVector(vIn);
                FIrrTimesVectorPerSubdomain.Add(FIrrTimesVectorOfSubdomain);
                FIrcTransposedTimesVectorPerSubdomain.Add(FIrcTransposedTimesVectorOfSubdomain);
            }

            Vector FIrrTimesVectorTotal = transferrer.SumVectors(FIrrTimesVectorPerSubdomain);
            Vector FIrcTransposedTimesVectorTotal = transferrer.SumVectors(FIrcTransposedTimesVectorPerSubdomain);
            return (FIrrTimesVectorTotal, FIrcTransposedTimesVectorTotal);
        }
    }
}