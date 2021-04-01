using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    internal class FetiDPGlobalDofSeparator
    {
        private readonly IModel model;

        internal FetiDPGlobalDofSeparator(IModel model)
        {
            this.model = model;
        }

        /// <summary>
        /// Dofs where Lagrange multipliers will be applied. These do not include corner dofs.
        /// </summary>
        internal Dictionary<INode, IDofType[]> GlobalBoundaryDofs { get; private set; }

        /// <summary>
        /// Dof ordering for corner dofs of the model: Each (INode, IDofType) pair is associated with the index of that dof into 
        /// a vector corresponding to all corner dofs of the model.
        /// </summary>
        internal DofTable GlobalCornerDofOrdering { get; private set; } // Only for master

        /// <summary>
        /// If Xf is a vector with all free dofs of the model and Xc is a vector with all corner dofs of the model, then
        /// Xf[GlobalCornerToFreeDofMap[i]] = Xc[i].
        /// </summary>
        internal int[] GlobalCornerToFreeDofMap { get; set; } // Only for master

        /// <summary>
        /// The number of corner dofs of the model.
        /// </summary>
        internal int NumGlobalCornerDofs { get; private set; }

        /// <summary>
        /// Bc unsigned boolean matrices that map global to subdomain corner dofs. This method must be called after 
        /// <see cref="DefineGlobalCornerDofs(Dictionary{int, HashSet{INode}})"/>.
        /// </summary>
        internal Dictionary<ISubdomain, UnsignedBooleanMatrix> CalcCornerMappingMatrices(
            Dictionary<ISubdomain, DofTable> subdomainCornerDofOrderings)
        { //TODO: Can I reuse subdomain data? Yes if the global corner dofs have not changed.
            var cornerBooleanMatrices =  new Dictionary<ISubdomain, UnsignedBooleanMatrix>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                DofTable localCornerDofOrdering = subdomainCornerDofOrderings[subdomain];
                int numLocalCornerDofs = localCornerDofOrdering.EntryCount;
                var Bc = new UnsignedBooleanMatrix(numLocalCornerDofs, NumGlobalCornerDofs);
                foreach ((INode node, IDofType dofType, int localIdx) in localCornerDofOrdering)
                {
                    int globalIdx = GlobalCornerDofOrdering[node, dofType];
                    Bc.AddEntry(localIdx, globalIdx);
                }
                cornerBooleanMatrices[subdomain] = Bc;
            }
            return cornerBooleanMatrices;
        }

        internal void DefineGlobalBoundaryDofs(HashSet<INode> globalCornerNodes)
        {
            IEnumerable<INode> globalRemainderNodes = model.EnumerateNodes().Where(node => !globalCornerNodes.Contains(node));
            GlobalBoundaryDofs =
                DofSeparationUtilities.DefineGlobalBoundaryDofs(globalRemainderNodes, model.GlobalDofOrdering.GlobalFreeDofs); //TODO: This could be reused in some cases
        }

        internal void DefineGlobalCornerDofs(HashSet<INode> globalCornerNodes)
        {
            // Order global corner dofs and create the global corner to global free map.
            var cornerToGlobalDofs = new List<int>(globalCornerNodes.Count * 3);
            var globalCornerDofOrdering = new DofTable(); //TODO: Should this be cached?
            int cornerDofCounter = 0;
            foreach (INode cornerNode in new SortedSet<INode>(globalCornerNodes)) //TODO: Must they be sorted?
            {
                bool hasFreeDofs = model.GlobalDofOrdering.GlobalFreeDofs.TryGetDataOfRow(cornerNode,
                    out IReadOnlyDictionary<IDofType, int> dofsOfNode);
                if (!hasFreeDofs) throw new Exception($"Corner node {cornerNode.ID} has only constrained or embedded dofs.");
                foreach (var dofTypeIdxPair in dofsOfNode)
                {
                    IDofType dofType = dofTypeIdxPair.Key;
                    int globalDofIdx = dofTypeIdxPair.Value;
                    globalCornerDofOrdering[cornerNode, dofType] = cornerDofCounter++;
                    cornerToGlobalDofs.Add(globalDofIdx);
                }
            }
            NumGlobalCornerDofs = cornerDofCounter;
            GlobalCornerToFreeDofMap = cornerToGlobalDofs.ToArray();
            GlobalCornerDofOrdering = globalCornerDofOrdering;
        }

        internal void ReorderGlobalCornerDofs(DofPermutation permutation)
        {
            if (permutation.IsBetter)
            {
                GlobalCornerDofOrdering.Reorder(permutation.PermutationArray, permutation.PermutationIsOldToNew);
                GlobalCornerToFreeDofMap = permutation.ReorderKeysOfDofIndicesMap(GlobalCornerToFreeDofMap);
            }
        }
    }
}
