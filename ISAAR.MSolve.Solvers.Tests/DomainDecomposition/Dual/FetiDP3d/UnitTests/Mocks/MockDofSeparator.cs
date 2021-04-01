using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks
{
    public class MockDofSeparator : IFetiDPDofSeparator
    {
        public DofTable GlobalCornerDofOrdering => throw new NotImplementedException();

        public int[] GlobalCornerToFreeDofMap => throw new NotImplementedException();

        public int NumGlobalCornerDofs => throw new NotImplementedException();

        public Dictionary<INode, IDofType[]> GlobalBoundaryDofs => throw new NotImplementedException();

        public int[] GetBoundaryDofIndices(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public (INode node, IDofType dofType)[] GetBoundaryDofs(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public UnsignedBooleanMatrix GetCornerBooleanMatrix(ISubdomain subdomain)
        {
            Matrix dense = ExpectedConnectivityData.GetMatrixBc(subdomain.ID);
            var sparse = new UnsignedBooleanMatrix(dense.NumRows, dense.NumColumns);
            for (int j = 0; j < dense.NumColumns; ++j)
            {
                for (int i = 0; i < dense.NumRows; ++i)
                {
                    if (dense[i, j] == 1) sparse.AddEntry(i, j);
                }
            }
            return sparse;
        }

        public int[] GetCornerDofIndices(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public DofTable GetCornerDofOrdering(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public IReadOnlyList<(INode node, IDofType dofType)> GetCornerDofs(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public int[] GetInternalDofIndices(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public int[] GetRemainderDofIndices(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public DofTable GetRemainderDofOrdering(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public void ReorderInternalDofs(IFetiDPSeparatedDofReordering reordering)
        {
            throw new NotImplementedException();
        }

        public void SeparateDofs(IFetiDPSeparatedDofReordering reordering)
        {
            throw new NotImplementedException();
        }
    }
}
