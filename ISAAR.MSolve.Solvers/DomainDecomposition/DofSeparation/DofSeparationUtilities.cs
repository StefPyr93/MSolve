using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation
{
    public static class DofSeparationUtilities
    {
        public static Dictionary<INode, IDofType[]> DefineGlobalBoundaryDofs(IEnumerable<INode> allNodes, 
            DofTable globalFreeDofs)
        {
            var globalBoundaryDofs = new Dictionary<INode, IDofType[]>();

            //TODO: model.Nodes probably doesn't work if there are embedded nodes. It is time to isolate the embedded nodes. Or I could use the GlobalDofOrdering.
            foreach (INode node in allNodes)
            {
                int nodeMultiplicity = node.Multiplicity;
                if (nodeMultiplicity > 1)
                {
                    // Access the free dofs only. Does this also filter out embedded dofs?
                    IDofType[] dofsOfNode = globalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.

                    // If all dofs of this node are constrained, then it is not considered boundary.
                    if (dofsOfNode.Length == 0) continue;
                    else globalBoundaryDofs[node] = dofsOfNode;
                }
            }

            return globalBoundaryDofs;
        }

        public static (int[] internalDofIndices, int[] boundaryDofIndices, (INode node, IDofType dofType)[] boundaryDofs)
            SeparateBoundaryInternalDofs(IEnumerable<INode> nodes, DofTable freeDofs)
        {
            var boundaryDofs = new SortedDictionary<int, (INode node, IDofType dofType)>(); // key = dofIdx, value = (node, dofType)
            var internalDofs = new SortedSet<int>(); //TODO: Why set instead of List?
            foreach (INode node in nodes)
            {
                int nodeMultiplicity = node.Multiplicity;
                if (nodeMultiplicity > 1) // boundary node
                {
                    bool isNodeFree = freeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<IDofType, int> dofTypesIndices); // This avoids embedded and constrained dofs
                    if (isNodeFree)
                    {
                        foreach (var dofTypeIdxPair in dofTypesIndices)
                        {
                            boundaryDofs.Add(dofTypeIdxPair.Value, (node, dofTypeIdxPair.Key));
                        }
                    }
                }
                else // internal nodes
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) internalDofs.Add(dof);
                }
            }

            // The following are sorted in increasing order of boundary dof indices
            return (internalDofs.ToArray(), boundaryDofs.Keys.ToArray(), boundaryDofs.Values.ToArray());
        }

        /// <summary>
        /// For debugging purposes.
        /// </summary>
        /// <param name="path"></param>
        /// <param name="separatedDofs"></param>
        public static void WriteSeparation(Dictionary<string, int[]> separatedDofs, string path, bool append)
        {
            using (var writer = new System.IO.StreamWriter(path, append))
            {
                foreach (string category in separatedDofs.Keys)
                {
                    writer.Write(category + ": ");
                    foreach (int dof in separatedDofs[category]) writer.Write(dof + " ");
                    writer.WriteLine();
                }
            }
        }
    }
}
