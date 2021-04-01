using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class LagrangeMultipliersEnumeratorMpi : ILagrangeMultipliersEnumerator
    {
        private readonly ICrosspointStrategy crosspointStrategy;
        private readonly IDofSeparator dofSeparator;
        private readonly LagrangeMultiplierSerializer lagrangeSerializer;
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        private List<LagrangeMultiplier> lagrangeMultipliers_master;
        private Dictionary<int, SignedBooleanMatrixColMajor> subdomainBooleanMatrices;

        public LagrangeMultipliersEnumeratorMpi(ProcessDistribution processDistribution, IModel model,
            ICrosspointStrategy crosspointStrategy, IDofSeparator dofSeparator)
        {
            this.procs = processDistribution;
            this.model = model;
            this.crosspointStrategy = crosspointStrategy;
            this.dofSeparator = dofSeparator;
            this.lagrangeSerializer = new LagrangeMultiplierSerializer(model.DofSerializer);
        }

        public IReadOnlyList<LagrangeMultiplier> LagrangeMultipliers
        {
            get
            {
                procs.CheckProcessIsMaster();
                return lagrangeMultipliers_master;
            }
        }

        public int NumLagrangeMultipliers { get; private set; }

        public SignedBooleanMatrixColMajor GetBooleanMatrix(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainBooleanMatrices[subdomain.ID];
        }

        public void CalcBooleanMatrices(Func<ISubdomain, DofTable> getSubdomainDofOrdering)
        {
            // Define the lagrange multipliers, serialize and broadcast them to other processes
            int[] serializedLagranges = null;
            if (procs.IsMasterProcess)
            {
                lagrangeMultipliers_master = LagrangeMultipliersUtilities.DefineLagrangeMultipliers(
                    dofSeparator.GlobalBoundaryDofs, crosspointStrategy);
                NumLagrangeMultipliers = lagrangeMultipliers_master.Count;
                serializedLagranges = lagrangeSerializer.Serialize(lagrangeMultipliers_master);
            }
            MpiUtilities.BroadcastArray(procs.Communicator, ref serializedLagranges, procs.MasterProcess);

            // Deserialize the lagrange multipliers in other processes and calculate the boolean matrices for each subdomain
            //ISubdomain subdomain = model.GetSubdomain(procs.OwnSubdomainID);
            //DofTable subdomainDofOrdering = getSubdomainDofOrdering(subdomain);
            var subdomainDofOrderings = new Dictionary<int, DofTable>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomainDofOrderings[s] = getSubdomainDofOrdering(subdomain);
            }
            if (procs.IsMasterProcess)
            {
                subdomainBooleanMatrices = CalcSubdomainBooleanMatrices(lagrangeMultipliers_master, subdomainDofOrderings);
            }
            else
            {
                //TODO: Using incomplete data and especially passing it between multiple classes is very fragile.
                Dictionary<int, ISubdomain> processSubdomains = procs.GetSubdomainsOfProcess(model);
                LagrangeMultiplier[] incompleteLagranges = lagrangeSerializer.DeserializeIncompletely(serializedLagranges,
                    processSubdomains);
                NumLagrangeMultipliers = incompleteLagranges.Length;
                subdomainBooleanMatrices = CalcSubdomainBooleanMatricesFromIncompleteData(incompleteLagranges, 
                    subdomainDofOrderings);
            }
        }

        // Optimized version of CalcSubdomainBooleanMatricesFromIncompleteData. Avoids null checks
        private static Dictionary<int, SignedBooleanMatrixColMajor> CalcSubdomainBooleanMatrices(
            IReadOnlyList<LagrangeMultiplier> globalLagranges, Dictionary<int, DofTable> remainderDofOrderings)
        {
            int numGlobalLagranges = globalLagranges.Count;
            var booleanMatrices = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (var subdomainOrderingPair in remainderDofOrderings)
            {
                int sub = subdomainOrderingPair.Key;
                DofTable ordering = subdomainOrderingPair.Value;
                booleanMatrices[sub] = new SignedBooleanMatrixColMajor(numGlobalLagranges, ordering.EntryCount);
            }

            for (int i = 0; i < numGlobalLagranges; ++i) // Global lagrange multiplier index
            {
                LagrangeMultiplier lagrange = globalLagranges[i];

                // Subdomain plus
                bool isInProcess = remainderDofOrderings.TryGetValue(lagrange.SubdomainPlus.ID, out DofTable orderingPlus);
                if (isInProcess)
                {
                    int dofIdx = orderingPlus[lagrange.Node, lagrange.DofType];
                    booleanMatrices[lagrange.SubdomainPlus.ID].AddEntry(i, dofIdx, true);
                }

                // Subdomain minus
                isInProcess = remainderDofOrderings.TryGetValue(lagrange.SubdomainMinus.ID, out DofTable orderingMinus);
                if (isInProcess)
                {
                    int dofIdx = orderingMinus[lagrange.Node, lagrange.DofType];
                    booleanMatrices[lagrange.SubdomainMinus.ID].AddEntry(i, dofIdx, false);
                }
            }

            return booleanMatrices;
        }

        private static Dictionary<int, SignedBooleanMatrixColMajor> CalcSubdomainBooleanMatricesFromIncompleteData(
            IReadOnlyList<LagrangeMultiplier> incompleteGlobalLagranges, Dictionary<int, DofTable> remainderDofOrderings)
        {
            int numGlobalLagranges = incompleteGlobalLagranges.Count;
            var booleanMatrices = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (var subdomainOrderingPair in remainderDofOrderings)
            {
                int sub = subdomainOrderingPair.Key;
                DofTable ordering = subdomainOrderingPair.Value;
                booleanMatrices[sub] = new SignedBooleanMatrixColMajor(numGlobalLagranges, ordering.EntryCount);
            }

            for (int i = 0; i < numGlobalLagranges; ++i) // Global lagrange multiplier index
            {
                LagrangeMultiplier lagrange = incompleteGlobalLagranges[i];
                if (lagrange == null) continue; // This lagrange's data is irrelevant and unavailable for the current subdomains.

                // Subdomain plus
                if (lagrange.SubdomainPlus != null)
                {
                    int sub = lagrange.SubdomainPlus.ID;
                    bool isInProcess = remainderDofOrderings.TryGetValue(sub, out DofTable ordering);
                    if (isInProcess)
                    {
                        int dofIdx = ordering[lagrange.Node, lagrange.DofType];
                        booleanMatrices[sub].AddEntry(i, dofIdx, true);
                    }
                }

                // Subdomain minus
                if (lagrange.SubdomainMinus != null)
                {
                    int sub = lagrange.SubdomainMinus.ID;
                    bool isInProcess = remainderDofOrderings.TryGetValue(sub, out DofTable ordering);
                    if (isInProcess)
                    {
                        int dofIdx = ordering[lagrange.Node, lagrange.DofType];
                        booleanMatrices[sub].AddEntry(i, dofIdx, false);
                    }
                }
            }

            return booleanMatrices;
        }
    }
}
