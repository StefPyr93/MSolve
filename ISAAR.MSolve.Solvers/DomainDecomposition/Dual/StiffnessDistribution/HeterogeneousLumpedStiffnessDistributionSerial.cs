using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HeterogeneousLumpedStiffnessDistributionSerial : IStiffnessDistribution
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiDPMatrixManager matrixManager;
        private readonly IHeterogeneousDistributionLoadScaling loadScaling;
        private readonly IModel model;

        private Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses;
        private Dictionary<ISubdomain, double[]> boundaryRelativeStiffnesses;
        private DiagonalMatrix Dlambda;

        public HeterogeneousLumpedStiffnessDistributionSerial(IModel model, IFetiDPDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPMatrixManager matrixManager,
            IHeterogeneousDistributionLoadScaling loadScaling)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.matrixManager = matrixManager;
            this.loadScaling = loadScaling;
            this.boundaryRelativeStiffnesses = new Dictionary<ISubdomain, double[]>();
        }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain) => boundaryRelativeStiffnesses[subdomain];

        public IMappingMatrix CalcBoundaryPreconditioningSignedBooleanMatrix(ILagrangeMultipliersEnumerator lagrangeEnumerator,
            ISubdomain subdomain, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
        {
            return new HeterogeneousStiffnessDistributionUtilities.ScalingBooleanMatrixImplicit(
                dofSeparator, subdomain, boundaryDofStiffnesses, Dlambda, boundarySignedBooleanMatrix);
        }

        public double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load) 
            => loadScaling.ScaleNodalLoad(subdomain, load, boundaryDofStiffnesses);

        public void Update()
        {
            //TODO: When can I reuse data?
            //if (subdomain.ConnectivityModified) //TODO: Is this what I should check?
            //{ }

            var matricesKff = new Dictionary<ISubdomain, IIndexable2D>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                matricesKff[subdomain] = matrixManager.GetSubdomainMatrixManager(subdomain).LinearSystem.Matrix;
            }

            this.boundaryDofStiffnesses =
                HeterogeneousStiffnessDistributionUtilities.CalcBoundaryDofStiffnesses(dofSeparator, matricesKff);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                boundaryRelativeStiffnesses[subdomain] = HeterogeneousStiffnessDistributionUtilities.CalcBoundaryDofCoefficients(
                    dofSeparator, subdomain, boundaryDofStiffnesses);
            }
            this.Dlambda = HeterogeneousStiffnessDistributionUtilities.BuildDlambda(
                lagrangesEnumerator, boundaryDofStiffnesses);
        }
    }
}
