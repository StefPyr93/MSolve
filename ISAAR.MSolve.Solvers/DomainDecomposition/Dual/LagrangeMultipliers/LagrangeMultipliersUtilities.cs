using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;

//TODO: Rename the "remainder" in most method arguments. If anything it should be "free/subdomain". Also avoid passing the entry count.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    /// <summary>
    /// Many of these methods overlap, because they have optimizations for serial / MPI environments.
    /// </summary>
    public static class LagrangeMultipliersUtilities
    {
        //TODO: If LagrangeMultiplier stored ids of node, dofType, subdomains then this methods would not be needed. 
        //      The overload that uses IReadOnlyList<LagrangeMultiplier> could be used instead. But that would require accessing
        //      INode and IDofType from ISubdomain's dictionaries (they are wrapped by accessors) which is much slower. The extra 
        //      cost would not matter in MPI implementations, since that accesses would be also done in processes other than 
        //      master to create the list of SubdomainLagrangeMultiplier after broadcasting the lagrangeMultipliers.
        //      Another problem is that IModel and ISubdomain would have to be updated to treat IDofType as an entity.
        //TODO: Another option that avoids messing with LagrangeMultiplier and IModel, ISubdomain is to never have an explicit
        //      IReadOnlyList<LagrangeMultiplier> in MPI implementation. That data would be stored as an int[] with each lagrange
        //      multiplier occupying 4 consecutive entries. Reading/writing the entries would be delegated to ILagrangeSerializer.
        //      Then I guess, this overload would be the only one actually used.
        public static SignedBooleanMatrixColMajor CalcSubdomainBooleanMatrixOLD(int numGlobalLagranges,
            IEnumerable<SubdomainLagrangeMultiplier> subdomainLagranges, DofTable subdomainDofOrdering)
        {
            int numSubdomainDofs = subdomainDofOrdering.EntryCount;
            var booleanMatrix = new SignedBooleanMatrixColMajor(numGlobalLagranges, numSubdomainDofs);
            foreach (SubdomainLagrangeMultiplier lagrange in subdomainLagranges)
            {
                int dofIdx = subdomainDofOrdering[lagrange.Node, lagrange.DofType];
                booleanMatrix.AddEntry(lagrange.GlobalLagrangeIndex, dofIdx, lagrange.SubdomainSign);
            }
            return booleanMatrix;
        }

        public static SignedBooleanMatrixColMajor CalcSubdomainBooleanMatrixOLD(ISubdomain subdomain, 
            IReadOnlyList<LagrangeMultiplier> globalLagranges, DofTable subdomainDofOrdering)
        {
            int numGlobalLagranges = globalLagranges.Count;
            int numSubdomainDofs = subdomainDofOrdering.EntryCount;
            var booleanMatrix = new SignedBooleanMatrixColMajor(numGlobalLagranges, numSubdomainDofs);

            for (int i = 0; i < numGlobalLagranges; ++i) // Global lagrange multiplier index
            {
                LagrangeMultiplier lagrange = globalLagranges[i];
                if (lagrange.SubdomainPlus.ID == subdomain.ID)
                {
                    int dofIdx = subdomainDofOrdering[lagrange.Node, lagrange.DofType];
                    booleanMatrix.AddEntry(i, dofIdx, true);
                }
                else if (lagrange.SubdomainMinus.ID == subdomain.ID)
                {
                    int dofIdx = subdomainDofOrdering[lagrange.Node, lagrange.DofType];
                    booleanMatrix.AddEntry(i, dofIdx, false);
                }
            }

            return booleanMatrix;
        }

        public static (Dictionary<int, SignedBooleanMatrixColMajor> booleanMatrices, LagrangeMultiplier[] lagranges)
            CalcBooleanMatricesAndLagrangesOLD(IModel model, int numLagrangeMultipliers,
            List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)> boundaryNodeData,
            Dictionary<int, int> numRemainderDofs, Dictionary<int, DofTable> remainderDofOrderings) //TODO: Rename the "remainder"
        {
            // Initialize the signed boolean matrices.
            var booleanMatrices = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                booleanMatrices[subdomain.ID] =
                    new SignedBooleanMatrixColMajor(numLagrangeMultipliers, numRemainderDofs[subdomain.ID]);
            }
            var lagrangeMultipliers = new LagrangeMultiplier[numLagrangeMultipliers];

            // Fill the boolean matrices and lagrange multiplier data: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int lag = 0; // Lagrange multiplier index
            foreach ((INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)
                in boundaryNodeData)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrixColMajor booleanPlus = booleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrixColMajor booleanMinus = booleanMatrices[subdomainsMinus[c].ID];

                    //TODO: The dof indices have already been accessed. Reuse it if possible.
                    IReadOnlyDictionary<IDofType, int> dofsPlus = remainderDofOrderings[subdomainsPlus[c].ID].GetDataOfRow(node);
                    IReadOnlyDictionary<IDofType, int> dofsMinus = remainderDofOrderings[subdomainsMinus[c].ID].GetDataOfRow(node);

                    foreach (IDofType dof in dofs)
                    {
                        booleanPlus.AddEntry(lag, dofsPlus[dof], true);
                        booleanMinus.AddEntry(lag, dofsMinus[dof], false);
                        lagrangeMultipliers[lag] = new LagrangeMultiplier(node, dof, subdomainsPlus[c], subdomainsMinus[c]);
                        ++lag;
                    }
                }
            }

            return (booleanMatrices, lagrangeMultipliers);
        }

        public static (List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)> combos,
            int numLagranges) DefineLagrangeCombinationsOLD(IDofSeparator dofSeparator, ICrosspointStrategy crosspointStrategy)
        {
            // Find boundary dual nodes and dofs
            var boundaryNodeData = new List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus,
                ISubdomain[] subdomainsMinus)>(dofSeparator.GlobalBoundaryDofs.Count);

            // Find continuity equations.
            int numLagrangeMultipliers = 0;
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                IDofType[] dofsOfNode = nodeDofsPair.Value;
                ISubdomain[] nodeSubdomains = node.SubdomainsDictionary.Values.ToArray();
                ISubdomain[] subdomainsPlus, subdomainsMinus;
                int multiplicity = nodeSubdomains.Length;
                if (multiplicity == 2)
                {
                    subdomainsPlus = new ISubdomain[] { nodeSubdomains[0] };
                    subdomainsMinus = new ISubdomain[] { nodeSubdomains[1] };
                }
                else (subdomainsPlus, subdomainsMinus) = crosspointStrategy.FindSubdomainCombinations(nodeSubdomains);

                boundaryNodeData.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
                numLagrangeMultipliers += dofsOfNode.Length * subdomainsPlus.Length;
            }

            return (boundaryNodeData, numLagrangeMultipliers);
        }

        public static List<LagrangeMultiplier> DefineLagrangeMultipliers(Dictionary<INode, IDofType[]> globalBoundaryDofs,
            ICrosspointStrategy crosspointStrategy)
        {
            var lagranges = new List<LagrangeMultiplier>();
            foreach (KeyValuePair<INode, IDofType[]> nodeDofsPair in globalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                IDofType[] dofsOfNode = nodeDofsPair.Value;
                ISubdomain[] nodeSubdomains = node.SubdomainsDictionary.Values.ToArray();

                // Find the how many lagrange multipliers are needed for each dof of this node 
                // and between which subdomains they should be applied to 
                ISubdomain[] subdomainsPlus, subdomainsMinus;
                int multiplicity = nodeSubdomains.Length;
                if (multiplicity == 2)
                {
                    subdomainsPlus = new ISubdomain[] { nodeSubdomains[0] };
                    subdomainsMinus = new ISubdomain[] { nodeSubdomains[1] };
                }
                else (subdomainsPlus, subdomainsMinus) = crosspointStrategy.FindSubdomainCombinations(nodeSubdomains);

                // Add the lagrange multipliers of this node to the global list. The order is important: 
                // node major - subdomain combination medium - dof minor.
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    foreach (IDofType dof in dofsOfNode)
                    {
                        lagranges.Add(new LagrangeMultiplier(node, dof, subdomainsPlus[c], subdomainsMinus[c]));
                    }
                }
            }
            return lagranges;
        }
    }
}
