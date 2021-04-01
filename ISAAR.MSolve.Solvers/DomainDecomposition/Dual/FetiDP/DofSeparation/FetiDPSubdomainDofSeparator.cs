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
    //TODO: This should be decoupled from MPI logic, so that it can be used in serial code too.
    internal class FetiDPSubdomainDofSeparator
    {
        internal FetiDPSubdomainDofSeparator(ISubdomain subdomain)
        {
            this.Subdomain = subdomain;
        }

        /// <summary>
        /// Indices of boundary remainder dofs into the sequence of all remainder dofs of a subdomain.
        /// </summary>
        internal int[] BoundaryDofIndices { get; private set; } //These can live inside each process. They are only needed globally if GatherGlobalDisplacements() is called (which normally won't be). In this case, they must be gathered.

        /// <summary>
        /// (Node, IDofType) pairs for each boundary remainder dof of a subdomain. Their order is the same one as 
        /// <see cref="BoundaryDofIndices"/>.
        /// </summary>
        internal (INode node, IDofType dofType)[] BoundaryDofs { get; private set; } //These are stored both in the process of each subdomain and in the master process.

        internal UnsignedBooleanMatrix CornerBooleanMatrix { get; private set; }

        /// <summary>
        /// Indices of (boundary) corner dofs into the sequence of all free dofs of a subdomain.
        /// </summary>
        internal int[] CornerDofIndices { get; private set; }

        /// <summary>
        /// Dof ordering for corner dofs of a subdomain: Each (INode, IDofType) pair of the subdomain is associated with the   
        /// index of that dof into a vector corresponding to corner dofs of that subdomain.
        /// </summary>
        internal DofTable CornerDofOrdering { get; private set; }

        /// <summary>
        /// Indices of internal remainder dofs into the sequence of all remainder dofs of a subdomain.
        /// </summary>
        internal int[] InternalDofIndices { get; private set; } 

        /// <summary>
        /// Dof ordering for remainder (boundary and internal) dofs of a subdomain: Each (INode, IDofType) pair of the 
        /// subdomain is associated with the index of that dof into a vector corresponding to remainder dofs of that subdomain.
        /// </summary>
        internal DofTable RemainderDofOrdering { get; private set; }

        /// <summary>
        /// Indices of remainder (boundary and internal) dofs into the sequence of all free dofs of a subdomain.
        /// </summary>
        internal int[] RemainderDofIndices { get; private set; }

        internal ISubdomain Subdomain { get; }

        internal IReadOnlyList<(INode node, IDofType dofType)> GetCornerDofs(HashSet<INode> cornerNodes)
        {
            // Optimization if the required capacity is known beforehand
            List<(INode node, IDofType dofType)> cornerDofs = null;
            if (CornerDofIndices != null) cornerDofs = new List<(INode node, IDofType dofType)>(CornerDofIndices.Length);
            else cornerDofs = new List<(INode node, IDofType dofType)>();

            foreach (INode node in cornerNodes)
            {
                foreach (IDofType dof in Subdomain.FreeDofOrdering.FreeDofs.GetColumnsOfRow(node)) cornerDofs.Add((node, dof));
            }

            return cornerDofs;
        }

        internal void ReorderInternalDofs(DofPermutation permutation)
        {
            if (permutation.IsBetter)
            {
                InternalDofIndices = permutation.ReorderKeysOfDofIndicesMap(InternalDofIndices);
            }
        }

        internal void ReorderRemainderDofs(DofPermutation permutation)
        {
            if (permutation.IsBetter)
            {
                RemainderDofIndices = permutation.ReorderKeysOfDofIndicesMap(RemainderDofIndices);
                RemainderDofOrdering.Reorder(permutation.PermutationArray, permutation.PermutationIsOldToNew);
            }
        }

        internal void SeparateBoundaryInternalDofs(HashSet<INode> cornerNodes)
        {
            IEnumerable<INode> remainderAndConstrainedNodes = 
                Subdomain.EnumerateNodes().Where(node => !cornerNodes.Contains(node));

            (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofConnectivities)
                = DofSeparationUtilities.SeparateBoundaryInternalDofs(remainderAndConstrainedNodes, RemainderDofOrdering);
            InternalDofIndices = internalDofIndices;
            BoundaryDofIndices = boundaryDofIndices;
            BoundaryDofs = boundaryDofConnectivities;
        }

        internal void SeparateCornerRemainderDofs(HashSet<INode> cornerNodes)
        {
            IEnumerable<INode> remainderAndConstrainedNodes = 
                Subdomain.EnumerateNodes().Where(node => !cornerNodes.Contains(node));

            var cornerDofs = new List<int>();
            var remainderDofs = new List<int>();
            foreach (INode node in cornerNodes)
            {
                IEnumerable<int> dofsOfNode = Subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                cornerDofs.AddRange(dofsOfNode);
            }
            foreach (INode node in remainderAndConstrainedNodes)
            {
                IEnumerable<int> dofsOfNode = Subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node);
                remainderDofs.AddRange(dofsOfNode);
            }
            CornerDofIndices = cornerDofs.ToArray();
            RemainderDofIndices = remainderDofs.ToArray();

            // This dof ordering will be optimized, such that the factorization of Krr is efficient.
            RemainderDofOrdering = Subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(remainderAndConstrainedNodes);
            CornerDofOrdering = Subdomain.FreeDofOrdering.FreeDofs.GetSubtableForNodes(cornerNodes);
        }

        internal void SetCornerBooleanMatrix(UnsignedBooleanMatrix matrix, object caller)
        {
            if (caller is IFetiDPDofSeparator) CornerBooleanMatrix = matrix;
            else throw new AccessViolationException("Caller cannot modify that memory");
        }
    }
}
