using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities
{
    public enum MatrixFormat
    {
        Dense, Skyline, SuiteSparse
    }

    public static class MatrixFormatSelection
    {
        public static IFetiDPMatrixManagerFactory DefineMatrixManagerFactory(MatrixFormat format)
        {
            if (format == MatrixFormat.Dense) return new FetiDPMatrixManagerFactoryDense();
            else if (format == MatrixFormat.Skyline) return new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
            else if (format == MatrixFormat.SuiteSparse) return new FetiDPMatrixManagerFactorySuitesparse();
            else throw new NotImplementedException();
        }
    }
}
