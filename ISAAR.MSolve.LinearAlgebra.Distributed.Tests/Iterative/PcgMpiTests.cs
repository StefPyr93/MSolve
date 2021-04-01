using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;
using MPI;
using ISAAR.MSolve.LinearAlgebra.Distributed.Iterative;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Iterative
{
    /// <summary>
    /// Tests for <see cref="PcgAlgorithm"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PcgMpiTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-5);

        public static void TestPosDefSystemOnMaster(int numProcesses)
        {
            Intracommunicator comm = Communicator.world;
            int master = 0;

            var builder = new PcgAlgorithmMpi.Builder();
            builder.ResidualTolerance = 1E-7;
            builder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1.0);
            var pcg = builder.Build(comm, master);

            Matrix A = null;
            JacobiPreconditioner M = null;
            if (comm.Rank == master)
            {
                A = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                M = new JacobiPreconditioner(A.GetDiagonalAsArray());
            }

            var matrix = new MasterOnlyMatrix(comm, master, SymmPosDef10by10.Order, A);
            var preconditioner = new MasterOnlyPreconditioner(comm, master, M);

            if (comm.Rank == master)
            {
                var b = Vector.CreateFromArray(SymmPosDef10by10.Rhs);
                Vector xComputed = Vector.CreateZero(A.NumRows);
                IterativeStatistics stats = pcg.Solve(matrix, preconditioner, b, xComputed, true, 
                    () => Vector.CreateZero(b.Length));

                var xExpected = Vector.CreateFromArray(SymmPosDef10by10.Lhs);
                comparer.AssertEqual(xExpected, xComputed);
            }
            else
            {
                IterativeStatistics stats = pcg.Solve(matrix, preconditioner, null, null, true, () => null);
            }
        }

        private class MasterOnlyMatrix : ILinearTransformation
        {
            private readonly Intracommunicator comm; //TODO: Use a dedicated class to handle this and master.
            private readonly int master;
            private readonly IMatrixView matrix_master;

            internal MasterOnlyMatrix(Intracommunicator comm, int masterProcess, int order, IMatrixView matrix_master)
            {
                this.comm = comm;
                this.master = masterProcess;
                this.matrix_master = matrix_master;

                this.NumRows = order;
                this.NumColumns = order;
            }

            public int NumColumns { get; }

            public int NumRows { get; }

            public void Multiply(IVectorView lhsVector, IVector rhsVector)
            {
                if (comm.Rank == master) matrix_master.MultiplyIntoResult(lhsVector, rhsVector);
            }
        }

        private class MasterOnlyPreconditioner : IPreconditioner
        {
            private readonly Intracommunicator comm; //TODO: Use a dedicated class to handle this and master.
            private readonly int master;
            private readonly IPreconditioner preconditioner_master;

            internal MasterOnlyPreconditioner(Intracommunicator comm, int masterProcess, IPreconditioner preconditioner_master)
            {
                this.comm = comm;
                this.master = masterProcess;
                this.preconditioner_master = preconditioner_master;

            }

            public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
            {
                if (comm.Rank == master) preconditioner_master.SolveLinearSystem(rhsVector, lhsVector);
            }
        }
    }
}
