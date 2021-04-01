using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

// 
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HomogeneousStiffnessDistributionSerial : IStiffnessDistribution
    {
        private readonly IFetiDPDofSeparator dofSeparator;

        /// <summary>
        /// Each process stores only the ones corresponding to its subdomain. Master stores all of them.
        /// </summary>
        private readonly Dictionary<ISubdomain, double[]> inverseBoundaryDofMultiplicities;

        private readonly IHomogeneousDistributionLoadScaling loadScaling;
        private readonly IModel model;

        public HomogeneousStiffnessDistributionSerial(IModel model, IFetiDPDofSeparator dofSeparator, 
            IHomogeneousDistributionLoadScaling loadScaling)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.loadScaling = loadScaling;
            this.inverseBoundaryDofMultiplicities = new Dictionary<ISubdomain, double[]>();
        }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain) => inverseBoundaryDofMultiplicities[subdomain];

        public IMappingMatrix CalcBoundaryPreconditioningSignedBooleanMatrix(ILagrangeMultipliersEnumerator lagrangeEnumerator,
            ISubdomain subdomain, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
        {
            return new HomogeneousStiffnessDistributionUtilities.ScalingBooleanMatrixImplicit(
                inverseBoundaryDofMultiplicities[subdomain], boundarySignedBooleanMatrix);
        }

        public double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load) => loadScaling.ScaleNodalLoad(subdomain, load);

        public void Update()
        {
            // Calculate and store the inverse boundary dof multiplicities of each subdomain
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.ConnectivityModified) //TODO: Is this what I should check?
                {
                    inverseBoundaryDofMultiplicities[subdomain] =
                        HomogeneousStiffnessDistributionUtilities.CalcBoundaryDofInverseMultiplicities(
                            subdomain, dofSeparator.GetBoundaryDofs(subdomain));
                }
            }
        }
    }
}