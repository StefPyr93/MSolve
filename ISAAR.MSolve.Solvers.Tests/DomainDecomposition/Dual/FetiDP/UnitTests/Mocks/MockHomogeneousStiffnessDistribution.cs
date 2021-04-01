using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests.Mocks
{
    public class MockHomogeneousStiffnessDistribution : IStiffnessDistribution
    {
        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain)
            => Example4x4QuadsHomogeneous.GetBoundaryDofCoefficients(subdomain.ID);

        public IMappingMatrix CalcBoundaryPreconditioningSignedBooleanMatrix(ILagrangeMultipliersEnumerator lagrangeEnumerator,
            ISubdomain subdomain, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
            => new MatrixBpbr(subdomain.ID);

        public double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load) => load.Amount / load.Node.Multiplicity;

        public void Update() { }

        private class MatrixBpbr : IMappingMatrix
        {
            private readonly Matrix Bpbr;

            internal MatrixBpbr(int subdomainID)
            {
                Bpbr = Example4x4QuadsHomogeneous.GetMatrixBpbr(subdomainID);
            }

            public int NumColumns => Bpbr.NumColumns;
            public int NumRows => Bpbr.NumRows;

            public Matrix CopyToFullMatrix() => Bpbr;

            public Vector Multiply(Vector vector, bool transposeThis = false) => Bpbr.Multiply(vector, transposeThis);

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false) => Bpbr.MultiplyRight(other, transposeThis);
        }
    }
}
