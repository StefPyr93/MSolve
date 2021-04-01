using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks
{
    public class MockMatrixManager : IFetiDP3dMatrixManager
    {
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager> subdomainMatrices;

        public MockMatrixManager(IModel model)
        {
            subdomainMatrices = new Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                subdomainMatrices[sub] = new MockSubdomainMatrixManager(sub);
            }
        }

        public Vector CoarseProblemRhs => ExpectedGlobalMatrices.VectorGlobalFcStarTildeSimple;

        public Vector GlobalDr => ExpectedGlobalMatrices.VectorDr;

        public void CalcCoarseProblemRhs() { }

        public void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection) { }

        public void ClearCoarseProblemRhs() { }

        public void ClearInverseCoarseProblemMatrix() { }

        public IFetiDPSubdomainMatrixManager GetFetiDPSubdomainMatrixManager(ISubdomain subdomain)
            => subdomainMatrices[subdomain];

        public IFetiSubdomainMatrixManager GetSubdomainMatrixManager(ISubdomain subdomain)
            => subdomainMatrices[subdomain];

        public Vector MultiplyInverseCoarseProblemMatrix(Vector vector)
        {
            if (vector != null) return ExpectedGlobalMatrices.MatrixGlobalKccStarTildeSimple.Invert() * vector;
            else return null;
        }

        public DofPermutation ReorderGlobalCornerDofs() => DofPermutation.CreateNoPermutation();

        public DofPermutation ReorderSubdomainInternalDofs(ISubdomain subdomain) => DofPermutation.CreateNoPermutation();

        public DofPermutation ReorderSubdomainRemainderDofs(ISubdomain subdomain) => DofPermutation.CreateNoPermutation();

        private class MockSubdomainMatrixManager : IFetiDPSubdomainMatrixManager
        {
            private readonly DiagonalMatrix invDii = null;
            private readonly Matrix invKii, invKrr, Kbb, Kbi, Kcc, Kff, Krc, Krr;

            internal MockSubdomainMatrixManager(ISubdomain subdomain)
            {
                int s = subdomain.ID;
                LinearSystem = new SingleSubdomainSystemMpi<Matrix>(subdomain);

                Fbc = ExpectedSubdomainMatrices.GetVectorFbc(s);
                FcStar = ExpectedSubdomainMatrices.GetVectorFcStar(s);
                Fr = ExpectedSubdomainMatrices.GetVectorFr(s);

                invDii = DiagonalMatrix.CreateFromArray(ExpectedSubdomainMatrices.GetMatrixKii(s).GetDiagonalAsArray());
                invDii.Invert();

                invKii = ExpectedSubdomainMatrices.GetMatrixKii(s).Invert();
                Kbb = ExpectedSubdomainMatrices.GetMatrixKbb(s);
                Kbi = ExpectedSubdomainMatrices.GetMatrixKbi(s);
                Kcc = ExpectedSubdomainMatrices.GetMatrixKcc(s);
                CoarseProblemSubmatrix = ExpectedSubdomainMatrices.GetMatrixKccStar(s);
                Kff = ExpectedSubdomainMatrices.GetMatrixKff(s);
                Krc = ExpectedSubdomainMatrices.GetMatrixKrc(s);
                Krr = ExpectedSubdomainMatrices.GetMatrixKrr(s);
                invKrr = ExpectedSubdomainMatrices.GetMatrixKrr(s).Invert();
            }

            public Vector Fbc { get; }
            public Vector FcStar { get; }
            public Vector Fr { get; }

            public IMatrixView CoarseProblemSubmatrix { get; }

            public ISingleSubdomainLinearSystemMpi LinearSystem { get; }

            public (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildFreeConstrainedMatrices(
                ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering, 
                IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
            {
                throw new NotImplementedException();
            }

            public void BuildFreeDofsMatrix(ISubdomainFreeDofOrdering dofOrdering, IElementMatrixProvider matrixProvider) { }

            public void CalcInverseKii(bool diagonalOnly) { }
            public void CalcSubdomainFcStartVector() { }
            public void CalcSubdomainKStarMatrices() { }

            public void ClearMatrices() { }
            public void ClearRhsVectors() { }
            public void CalcCoarseProblemSubmatrices() { }
            public void CalcCoarseProblemRhsSubvectors() { }

            public void ExtractBoundaryInternalSubmatricesAndInvertKii(bool diagonalKii) { }

            public void ExtractCornerRemainderRhsSubvectors() { }
            public void ExtractCornerRemainderSubmatrices() { }
            public void ExtractKbb() { }
            public void ExtractKbiKib() { }
            public void HandleDofOrderingWillBeModified() { }
            public void InvertKrr(bool inPlace) { }

            public Vector MultiplyInverseKiiTimes(Vector vector, bool diagonalOnly)
            {
                if (diagonalOnly) return invDii * vector;
                else return invKii * vector;
            }

            public Matrix MultiplyInverseKiiTimes(Matrix matrix, bool diagonalOnly)
            {
                if (diagonalOnly) return invDii * matrix;
                else return invKii * matrix;
            }

            public Vector MultiplyInverseKrrTimes(Vector vector) => invKrr * vector;
            public Vector MultiplyKbbTimes(Vector vector) => Kbb * vector;
            public Matrix MultiplyKbbTimes(Matrix matrix) => Kbb * matrix;
            public Vector MultiplyKbiTimes(Vector vector) => Kbi * vector;
            public Matrix MultiplyKbiTimes(Matrix matrix) => Kbi * matrix;
            public Vector MultiplyKccTimes(Vector vector) => Kcc * vector;
            public Vector MultiplyKcrTimes(Vector vector) => Krc.Multiply(vector, true);
            public Vector MultiplyKibTimes(Vector vector) => Kbi.Multiply(vector, true);
            public Matrix MultiplyKibTimes(Matrix matrix) => Kbi.MultiplyRight(matrix, true);
            public Vector MultiplyKrcTimes(Vector vector) => Krc * vector;

            public DofPermutation ReorderInternalDofs() => DofPermutation.CreateNoPermutation();
            public DofPermutation ReorderRemainderDofs() => DofPermutation.CreateNoPermutation();

            public IMatrixView CalcMatrixSb()
            {
                throw new NotImplementedException();
            }
        }
    }
}
