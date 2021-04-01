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

//TODO: This creates larger coarse problems and more MPI communication per PCG iteration. However the PCG iterations are fewer!!!
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public class AugmentationConstraintsRedundant : IAugmentationConstraints
    {
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IModel model;

        public AugmentationConstraintsRedundant(IModel model, IMidsideNodesSelection midsideNodesSelection,
            ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            this.model = model;
            this.MidsideNodesSelection = midsideNodesSelection;
            this.lagrangesEnumerator = lagrangesEnumerator;
        }

        public IMappingMatrix MatrixGlobalQr { get; private set; }

        public IMidsideNodesSelection MidsideNodesSelection { get; }

        public int NumGlobalAugmentationConstraints { get; private set; }

        public DofTable GlobalAugmentationDofOrdering => throw new NotImplementedException();

        public void CalcAugmentationMappingMatrices()
        {
            Table<INode, IDofType, HashSet<int>> augmentationLagranges = FindAugmentationLagranges();
            NumGlobalAugmentationConstraints = 0;
            foreach ((INode node, IDofType dof, HashSet<int> val) in augmentationLagranges)
            {
                NumGlobalAugmentationConstraints += val.Count;
            }

            var Qr = new UnsignedBooleanMatrix(lagrangesEnumerator.NumLagrangeMultipliers, NumGlobalAugmentationConstraints);
            int col = 0;
            foreach (INode node in MidsideNodesSelection.MidsideNodesGlobal)
            {
                foreach (IDofType dof in MidsideNodesSelection.DofsPerNode)
                {
                    foreach (int idx in augmentationLagranges[node, dof])
                    {
                        Qr.AddEntry(idx, col);
                        ++col;
                    }
                }
            }

            MatrixGlobalQr = Qr;
        }

        public DofTable GetAugmentationDofOrdering(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public GlobalToLocalBooleanMatrix GetMatrixBa(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public IMappingMatrix GetMatrixR1(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public int GetNumAugmentationDofs(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        private Table<INode, IDofType, HashSet<int>> FindAugmentationLagranges()
        {
            var augmentationLagranges = new Table<INode, IDofType, HashSet<int>>();
            foreach (INode node in MidsideNodesSelection.MidsideNodesGlobal)
            {
                foreach (IDofType dof in MidsideNodesSelection.DofsPerNode) augmentationLagranges[node, dof] = new HashSet<int>();
            }

            var midsideNodes = new HashSet<INode>(MidsideNodesSelection.MidsideNodesGlobal); // for faster look-ups. TODO: Use the table for look-ups
            for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; ++i)
            {
                LagrangeMultiplier lagr = lagrangesEnumerator.LagrangeMultipliers[i];
                if (midsideNodes.Contains(lagr.Node))
                {
                    Debug.Assert(MidsideNodesSelection.DofsPerNode.Contains(lagr.DofType));
                    augmentationLagranges[lagr.Node, lagr.DofType].Add(i);
                }
            }
            return augmentationLagranges;
        }

        //private static Dictionary<INode, List<int>> FindAugmentationLagranges(IMidsideNodesSelection midsideNodesSelection,
        //    ILagrangeMultipliersEnumerator lagrangesEnumerator)
        //{
        //    var augmentationLagranges = new Dictionary<INode, List<int>>();
        //    foreach (INode node in midsideNodesSelection.MidsideNodesGlobal) augmentationLagranges[node] = new List<int>();

        //    for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; ++i)
        //    {
        //        LagrangeMultiplier lagr = lagrangesEnumerator.LagrangeMultipliers[i];
        //        if (augmentationLagranges.ContainsKey(lagr.Node))
        //        {
        //            if (lagr.Node.ID == 38)
        //            {
        //                Console.WriteLine();
        //            }

        //            augmentationLagranges[lagr.Node].Add(i);
        //        }
        //    }
        //    return augmentationLagranges;
        //}

        public class Factory : IAugmentationConstraintsFactory
        {
            public IAugmentationConstraints CreateAugmentationConstraints(IModel model, 
                IMidsideNodesSelection midsideNodesSelection, IFetiDPDofSeparator dofSeparator,
                ILagrangeMultipliersEnumerator lagrangesEnumerator)
            {
                return new AugmentationConstraintsRedundant(model, midsideNodesSelection, lagrangesEnumerator);                               
            }
        }

    }
}
