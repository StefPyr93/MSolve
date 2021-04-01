using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public class AugmentationConstraints : IAugmentationConstraints
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IModel model;

        private Dictionary<ISubdomain, GlobalToLocalBooleanMatrix> matricesBa = 
            new Dictionary<ISubdomain, GlobalToLocalBooleanMatrix>();

        private Dictionary<ISubdomain, LocalToGlobalMappingMatrix> matricesR1 = 
            new Dictionary<ISubdomain, LocalToGlobalMappingMatrix>();

        private UnsignedBooleanMatrixColMajor matrixQr;

        private Dictionary<ISubdomain, DofTable> subdomainAugmentationDofOrderings =
            new Dictionary<ISubdomain, DofTable>();

        public AugmentationConstraints(IModel model, IMidsideNodesSelection midsideNodesSelection,
            IFetiDPDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            this.model = model;
            this.MidsideNodesSelection = midsideNodesSelection;
            this.dofSeparator = dofSeparator;
            this.lagrangesEnumerator = lagrangesEnumerator;
        }

        public DofTable GlobalAugmentationDofOrdering { get; private set; }

        public IMappingMatrix MatrixGlobalQr => matrixQr;

        public IMidsideNodesSelection MidsideNodesSelection { get; }

        public int NumGlobalAugmentationConstraints { get; private set; }

        public void CalcAugmentationMappingMatrices()
        {
            this.GlobalAugmentationDofOrdering = OrderGlobalAugmentationConstraints();
            CalcGlobalMatrixQr(GlobalAugmentationDofOrdering);
            CalcMatricesBa(GlobalAugmentationDofOrdering);
            CalcMatricesR1();
        }

        public DofTable GetAugmentationDofOrdering(ISubdomain subdomain) => subdomainAugmentationDofOrderings[subdomain];

        public GlobalToLocalBooleanMatrix GetMatrixBa(ISubdomain subdomain) => matricesBa[subdomain];

        public IMappingMatrix GetMatrixR1(ISubdomain subdomain) => matricesR1[subdomain];

        public int GetNumAugmentationDofs(ISubdomain subdomain) => matricesBa[subdomain].NumRows;

        private void CalcGlobalMatrixQr(DofTable augmentedOrdering)
        {
            matrixQr = 
                new UnsignedBooleanMatrixColMajor(lagrangesEnumerator.NumLagrangeMultipliers, NumGlobalAugmentationConstraints);
            for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; ++i)
            {
                LagrangeMultiplier lagr = lagrangesEnumerator.LagrangeMultipliers[i];
                bool isMidside = augmentedOrdering.TryGetValue(lagr.Node, lagr.DofType, out int augmentedIdx);
                if (isMidside) matrixQr.AddEntry(i, augmentedIdx);
            }
        }

        private void CalcMatricesBa(DofTable augmentedOrdering)
        {
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> midsideNodes = MidsideNodesSelection.GetMidsideNodesOfSubdomain(subdomain);
                int numSubdomainAugmentedConstraints = midsideNodes.Count * MidsideNodesSelection.DofsPerNode.Length;
                var augmentationDofOrdering = new DofTable();
                var subdomainToGlobalAugmentedConstraints = new int[numSubdomainAugmentedConstraints];

                int subdomainDofIdx = 0;
                foreach (INode node in midsideNodes)
                {
                    foreach (IDofType dof in MidsideNodesSelection.DofsPerNode)
                    {
                        int globalDofIdx = augmentedOrdering[node, dof];
                        augmentationDofOrdering[node, dof] = subdomainDofIdx;
                        subdomainToGlobalAugmentedConstraints[subdomainDofIdx] = globalDofIdx;
                        ++subdomainDofIdx;
                    }
                }

                this.subdomainAugmentationDofOrderings[subdomain] = augmentationDofOrdering;
                matricesBa[subdomain] =
                    new GlobalToLocalBooleanMatrix(NumGlobalAugmentationConstraints, subdomainToGlobalAugmentedConstraints);
            }
        }

        private void CalcMatricesR1()
        {
            // Q1[s] = Qr * Ba[s]^T
            // R1[s] = Br[s]^T * Q1[s] 
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                // Optimized sparse matrix multiplication for these 2 matrices (only in this case): The dot product of a column 
                // of Br and a column of Q1 will be nonzero, only if they both refer to the same midisde dof.
                SignedBooleanMatrixColMajor Br = lagrangesEnumerator.GetBooleanMatrix(subdomain);
                GlobalToLocalBooleanMatrix Ba = matricesBa[subdomain];
                
                // Find the columns of Br that correspond to the same dofs as columns of Q1. 
                // This must be done in the same order as they were created in Ba. 
                //TODO: Perhaps this should have bene stored somewehere.
                DofTable remainderDofs = dofSeparator.GetRemainderDofOrdering(subdomain);
                var midsideBrColumns = new int[Ba.NumRows];
                int localAugmentedDofIdx = 0;
                foreach (INode node in MidsideNodesSelection.GetMidsideNodesOfSubdomain(subdomain))
                {
                    foreach (IDofType dof in MidsideNodesSelection.DofsPerNode)
                    {
                        midsideBrColumns[localAugmentedDofIdx] = remainderDofs[node, dof];
                        ++localAugmentedDofIdx;
                    }
                }

                // Perform the matrix multiplication transpose(Br) * Q1: dot products of corresponding columns of Br and Q1
                // Columns of Q1 are read directly from Qr, by selecting the correct ones as specified by Ba
                int[] columnsOfQrToKeep = Ba.GetRowsToColumnsMap();
                var valuesR1 = new double[Ba.NumRows];
                for (int i = 0; i < Ba.NumRows; ++i)
                {
                    Dictionary<int, int> colBr = Br.GetColumn(midsideBrColumns[i]);
                    HashSet<int> colQ1 = matrixQr.GetNonZeroRowsOfColumn(columnsOfQrToKeep[i]);

                    // In general the column of Br has fewer (non-zero) entries than the column of Q1, thus iterating the former 
                    // could be slightly faster.
                    double sum = 0.0;
                    foreach (var rowSignPair in colBr)
                    {
                        if (colQ1.Contains(rowSignPair.Key)) sum += rowSignPair.Value;
                    }
                    valuesR1[i] = sum;
                }
                matricesR1[subdomain] = new LocalToGlobalMappingMatrix(Br.NumColumns, valuesR1, midsideBrColumns);
            }
        }

        private DofTable OrderGlobalAugmentationConstraints()
        {
            this.NumGlobalAugmentationConstraints =
                MidsideNodesSelection.DofsPerNode.Length * MidsideNodesSelection.MidsideNodesGlobal.Count;

            var ordering = new DofTable();
            int idx = 0;
            foreach (INode node in MidsideNodesSelection.MidsideNodesGlobal)
            {
                foreach (IDofType dof in MidsideNodesSelection.DofsPerNode) ordering[node, dof] = idx++;
            }
            return ordering;
        }

        public class Factory : IAugmentationConstraintsFactory
        {
            public IAugmentationConstraints CreateAugmentationConstraints(IModel model, 
                IMidsideNodesSelection midsideNodesSelection, IFetiDPDofSeparator dofSeparator, 
                ILagrangeMultipliersEnumerator lagrangesEnumerator)
            {
                return new AugmentationConstraints(model, midsideNodesSelection, dofSeparator, lagrangesEnumerator);
            }
        }
    }
}
