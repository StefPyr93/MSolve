using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Perhaps I should make this a static utility class
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public abstract class HomogeneousStiffnessDistributionOLD : IStiffnessDistributionOLD
    {
        private readonly IDofSeparator dofSeparator;
        private readonly Dictionary<int, double[]> inverseBoundaryDofMultiplicities;
        private readonly IModel model;

        public HomogeneousStiffnessDistributionOLD(IModel model, IDofSeparator dofSeparator)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.inverseBoundaryDofMultiplicities = new Dictionary<int, double[]>();
        }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain) => inverseBoundaryDofMultiplicities[subdomain.ID];

        //public Dictionary<int, double> CalcBoundaryDofCoefficientsOLD(INode node, IDofType dofType)
        //{
        //    var coeffs = new Dictionary<int, double>();
        //    double inverseMultiplicity = 1.0 / node.Multiplicity;
        //    foreach (int subdomainID in node.SubdomainsDictionary.Keys) coeffs[subdomainID] = inverseMultiplicity;
        //    return coeffs;
        //}

        public Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            ILagrangeMultipliersEnumeratorOLD lagrangeEnumerator,
            Dictionary<int, SignedBooleanMatrixColMajor> boundarySignedBooleanMatrices)
        {
            var matricesBpb = new Dictionary<int, IMappingMatrix>();
            foreach (int s in boundarySignedBooleanMatrices.Keys)
            {
                matricesBpb[s] = new HomogeneousStiffnessDistributionUtilities.ScalingBooleanMatrixImplicit(
                    inverseBoundaryDofMultiplicities[s], boundarySignedBooleanMatrices[s]);
            }
            return matricesBpb;
        }

        public void Update(Dictionary<int, IMatrixView> stiffnessesFreeFree)
        {
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;
                if (model.GetSubdomain(s).ConnectivityModified)
                {
                    Debug.WriteLine(
                        $"{this.GetType().Name}: Calculating the inverse multiplicities of the boundary dofs of subdomain {s}");
                    (INode node, IDofType dofType)[] boundaryDofs = dofSeparator.GetBoundaryDofs(subdomain);
                    var inverseMultiplicities = new double[boundaryDofs.Length];
                    for (int i = 0; i < boundaryDofs.Length; ++i)
                    {
                        inverseMultiplicities[i] = 1.0 / boundaryDofs[i].node.Multiplicity;
                    }
                    this.inverseBoundaryDofMultiplicities[s] = inverseMultiplicities;
                }
            }
        }

        public abstract double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load);
    }
}