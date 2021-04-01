using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class LagrangeMultipliersEnumeratorSerial : ILagrangeMultipliersEnumerator
    {
        private readonly ICrosspointStrategy crosspointStrategy;
        public readonly IDofSeparator dofSeparator;
        private readonly IModel model;

        private Dictionary<ISubdomain, SignedBooleanMatrixColMajor> subdomainBooleanMatrices;

        public LagrangeMultipliersEnumeratorSerial(IModel model, ICrosspointStrategy crosspointStrategy,
            IDofSeparator dofSeparator)
        {
            this.model = model;
            this.crosspointStrategy = crosspointStrategy;
            this.dofSeparator = dofSeparator;
        }

        public IReadOnlyList<LagrangeMultiplier> LagrangeMultipliers { get; private set; }

        public int NumLagrangeMultipliers { get; private set; }

        public SignedBooleanMatrixColMajor GetBooleanMatrix(ISubdomain subdomain) => subdomainBooleanMatrices[subdomain];

        public void CalcBooleanMatrices(Func<ISubdomain, DofTable> getSubdomainDofOrdering)
        {
            // Define the lagrange multipliers
            LagrangeMultipliers = 
                LagrangeMultipliersUtilities.DefineLagrangeMultipliers(dofSeparator.GlobalBoundaryDofs, crosspointStrategy);
            NumLagrangeMultipliers = LagrangeMultipliers.Count;

            // Define the subdomain dofs
            var subdomainDofOrderings = new Dictionary<ISubdomain, DofTable>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainDofOrderings[subdomain] = getSubdomainDofOrdering(subdomain);
            }

            // Calculate the boolean matrices
            subdomainBooleanMatrices = CalcAllBooleanMatrices(LagrangeMultipliers, subdomainDofOrderings);
        }

        /// <summary>
        /// This method is slower than <see cref="CalcBooleanMatricesAndLagranges(IModel, int, 
        /// List{(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)}, Dictionary{int, int}, 
        /// Dictionary{int, DofTable})"/>. It probably does not matter that much though.
        /// </summary>
        private static Dictionary<ISubdomain, SignedBooleanMatrixColMajor> CalcAllBooleanMatrices(
            IReadOnlyList<LagrangeMultiplier> globalLagranges, Dictionary<ISubdomain, DofTable> subdomainDofOrderings)
        {
            // Initialize the signed boolean matrices
            var booleanMatrices = new Dictionary<ISubdomain, SignedBooleanMatrixColMajor>();
            foreach (ISubdomain subdomain in subdomainDofOrderings.Keys)
            {
                booleanMatrices[subdomain] =
                    new SignedBooleanMatrixColMajor(globalLagranges.Count, subdomainDofOrderings[subdomain].EntryCount);
            }

            // Fill all boolean matrices simultaneously
            for (int i = 0; i < globalLagranges.Count; ++i) // Global lagrange multiplier index
            {
                LagrangeMultiplier lagrange = globalLagranges[i];

                int dofIdxPlus = subdomainDofOrderings[lagrange.SubdomainPlus][lagrange.Node, lagrange.DofType];
                booleanMatrices[lagrange.SubdomainPlus].AddEntry(i, dofIdxPlus, true);

                int dofIdxMinus = subdomainDofOrderings[lagrange.SubdomainMinus][lagrange.Node, lagrange.DofType];
                booleanMatrices[lagrange.SubdomainMinus].AddEntry(i, dofIdxMinus, false);
            }

            return booleanMatrices;
        }
    }
}
