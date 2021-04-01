using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

// 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HomogeneousStiffnessDistributionMpi : IStiffnessDistribution
    {
        private readonly IFetiDPDofSeparator dofSeparator;

        /// <summary>
        /// Each process stores only the ones corresponding to their subdomains. Master may store all of them.
        /// </summary>
        private readonly Dictionary<int, double[]> inverseBoundaryDofMultiplicities;

        private readonly IHomogeneousDistributionLoadScaling loadScaling;
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public HomogeneousStiffnessDistributionMpi(ProcessDistribution processDistribution, IModel model,
            IFetiDPDofSeparator dofSeparator, IHomogeneousDistributionLoadScaling loadScaling)
        {
            this.procs = processDistribution;
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.loadScaling = loadScaling;
            this.inverseBoundaryDofMultiplicities = new Dictionary<int, double[]>();
        }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return inverseBoundaryDofMultiplicities[subdomain.ID];
        }

        public IMappingMatrix CalcBoundaryPreconditioningSignedBooleanMatrix(ILagrangeMultipliersEnumerator lagrangeEnumerator,
            ISubdomain subdomain, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return new HomogeneousStiffnessDistributionUtilities.ScalingBooleanMatrixImplicit(
                inverseBoundaryDofMultiplicities[subdomain.ID], boundarySignedBooleanMatrix);
        }

        /// <summary>
        /// This is not necessary in a typical execution. It would be useful e.g. to create a global vector with all 
        /// displacements.
        /// </summary>
        public void GatherDataInMaster()
        {
            var activeSubdomains = new ActiveSubdomains(procs, s => model.GetSubdomain(s).ConnectivityModified);

            // Gather all boundary dof multiplicites in master. It is faster to send int[] than double[].
            var transferrer = new TransferrerPerSubdomain(procs);
            GetArrayLengthOfPackedData<double[]> getPackedDataLength = (s, arry) => arry.Length;
            PackSubdomainDataIntoArray<double[], int> packData =
                (s, inverse, direct, offsetDirect) => Invert(inverse, direct, offsetDirect);
            UnpackSubdomainDataFromArray<double[], int> unpackData =
                (s, direct, start, end) => Invert(direct, start, end);
            var allData = transferrer.GatherFromSomeSubdomainsPacked(
                this.inverseBoundaryDofMultiplicities, getPackedDataLength, packData, unpackData, activeSubdomains);
            if (procs.IsMasterProcess)
            {
                foreach (var pair in allData) this.inverseBoundaryDofMultiplicities[pair.Key] = pair.Value;
            }
        }

        public double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load) => loadScaling.ScaleNodalLoad(subdomain, load);

        public void Update() 
        {
            // Calculate and store the inverse boundary dof multiplicities of each modified subdomain in this process
            foreach (ISubdomain subdomain in procs.GetSubdomainsOfProcess(model).Values)
            {
                if (subdomain.ConnectivityModified) //TODO: Is this what I should check?
                {
                    inverseBoundaryDofMultiplicities[subdomain.ID] =
                        HomogeneousStiffnessDistributionUtilities.CalcBoundaryDofInverseMultiplicities(
                            subdomain, dofSeparator.GetBoundaryDofs(subdomain));
                }
            }
        }
        
        private static int[] Invert(double[] inverse, int[] direct, int offsetDirect)
        {
            for (int i = 0; i < inverse.Length; ++i) direct[offsetDirect + i] = (int)(Math.Round(1.0 / inverse[i]));
            return direct;
        }

        private static double[] Invert(int[] direct, int start, int end)
        {
            int length = end - start;
            var inverse = new double[length];
            for (int i = 0; i < length; ++i) inverse[i] = 1.0 / direct[start + i];
            return inverse;
        }
    }
}