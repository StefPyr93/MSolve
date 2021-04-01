using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public class AugmentationLagrangesMapping
    {

        private (Dictionary<ISubdomain, UnsignedBooleanMatrix>, Dictionary<ISubdomain, int[]>) CalcMappingsAugmentedToRemainder(IModel model,
            IFetiDPDofSeparator dofSeparator,
            IMidsideNodesSelection midsideNodesSelection, IDofType[] dofsPerNode,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, Dictionary<ISubdomain, DofTable> subdomainDofOrderings)
        {
            var mappings = new Dictionary<ISubdomain, UnsignedBooleanMatrix>();
            Dictionary<ISubdomain, int[]> subdQrtoGlobalQr = new Dictionary<ISubdomain, int[]>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int numRows = dofSeparator.GetRemainderDofIndices(subdomain).Length;
                int numCols = midsideNodesSelection.GetMidsideNodesOfSubdomain(subdomain).Count * dofsPerNode.Length;
                mappings[subdomain] = new UnsignedBooleanMatrix(numRows, numCols);
                subdQrtoGlobalQr[subdomain] = new int[numCols];
            }

            Table<INode, IDofType, HashSet<int>> augmentationLagranges =
                FindAugmentationLagranges(midsideNodesSelection, dofsPerNode, lagrangesEnumerator);
            int NumGlobalAugmentationConstraints = dofsPerNode.Length * midsideNodesSelection.MidsideNodesGlobal.Count;

            


            Dictionary<ISubdomain, int> subdOffsets = new Dictionary<ISubdomain, int>();
            for (int n = 0; n < midsideNodesSelection.MidsideNodesGlobal.Count; ++n)
            {
                INode node = midsideNodesSelection.MidsideNodesGlobal[n];
                int offset = n * dofsPerNode.Length;
                for (int j = 0; j < dofsPerNode.Length; ++j)
                {
                    HashSet<int> rowIndices = augmentationLagranges[node, dofsPerNode[j]];
                    //foreach (int i in rowIndices) MatrixQr[i, offset + j] = 1.0;

                    foreach (int i in rowIndices)
                    {
                        LagrangeMultiplier lagr = lagrangesEnumerator.LagrangeMultipliers[i];


                        int dofIdxPlus = subdomainDofOrderings[lagr.SubdomainPlus][lagr.Node, lagr.DofType];
                        mappings[lagr.SubdomainPlus].AddEntry(dofIdxPlus, subdOffsets[lagr.SubdomainPlus]);
                        subdQrtoGlobalQr[lagr.SubdomainPlus][subdOffsets[lagr.SubdomainPlus]] = i;
                        subdOffsets[lagr.SubdomainPlus]++;

                    }

                }
            }


            return (mappings, subdQrtoGlobalQr);
        }

        private static Table<INode, IDofType, HashSet<int>> FindAugmentationLagranges(
            IMidsideNodesSelection midsideNodesSelection, IEnumerable<IDofType> dofsPerNode,
            ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            var augmentationLagranges = new Table<INode, IDofType, HashSet<int>>();
            foreach (INode node in midsideNodesSelection.MidsideNodesGlobal)
            {
                foreach (IDofType dof in dofsPerNode) augmentationLagranges[node, dof] = new HashSet<int>();
            }

            var midsideNodes = new HashSet<INode>(midsideNodesSelection.MidsideNodesGlobal); // for faster look-ups. TODO: Use the table for look-ups
            for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; ++i)
            {
                LagrangeMultiplier lagr = lagrangesEnumerator.LagrangeMultipliers[i];
                if (midsideNodes.Contains(lagr.Node))
                {
                    Debug.Assert(dofsPerNode.Contains(lagr.DofType));
                    augmentationLagranges[lagr.Node, lagr.DofType].Add(i);
                }
            }
            return augmentationLagranges;
        }
    }
}
