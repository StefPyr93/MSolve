using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Direct
{
	public class CSparseLUSolver : SingleSubdomainSolverBase<CscMatrix>
	{
		private readonly double factorizationPivotTolerance;

		private bool factorizeInPlace = true;
		private bool mustFactorize = true;
		private LUCSparseNet factorization;

		private CSparseLUSolver(IModel model, double factorizationPivotTolerance, IDofOrderer dofOrderer) :
			base(model, dofOrderer, new CscAssembler(), "CSparseLUSolver")
		{
			this.factorizationPivotTolerance = factorizationPivotTolerance;
		}

		public override void HandleMatrixWillBeSet()
		{
			mustFactorize = true;
			factorization = null;
		}

		public override void Initialize() { }

		public override void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

		/// <summary>
		/// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
		/// </summary>
		public override void Solve()
		{
			var watch = new Stopwatch();
			if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
			//else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

			// Factorization
			if (mustFactorize)
			{
				watch.Start();
				factorization = LUCSparseNet.Factorize(linearSystem.Matrix, factorizationPivotTolerance);
				watch.Stop();
				Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
				watch.Reset();
				mustFactorize = false;
			}

			// Substitutions
			watch.Start();
			factorization.SolveLinearSystem(linearSystem.RhsConcrete, linearSystem.SolutionConcrete);
			watch.Stop();
			Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
			Logger.IncrementAnalysisStep();

			Debug.WriteLine($"||res||={linearSystem.Matrix.Multiply(linearSystem.SolutionConcrete).Subtract(linearSystem.RhsConcrete).Norm2()}");
		}

		protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
		{
			var watch = new Stopwatch();

			// Factorization
			if (mustFactorize)
			{
				watch.Start();
				factorization = LUCSparseNet.Factorize(linearSystem.Matrix, factorizationPivotTolerance);
				watch.Stop();
				Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
				watch.Reset();
				mustFactorize = false;
			}

			// Substitutions
			watch.Start();

			// Solve each linear system separately, to avoid copying the RHS matrix to a dense one.
			int systemOrder = linearSystem.Matrix.NumColumns;
			int numRhs = otherMatrix.NumColumns;
			var solutionVectors = Matrix.CreateZero(systemOrder, numRhs);
			Vector solutionVector = linearSystem.CreateZeroVectorConcrete();
			for (int j = 0; j < numRhs; ++j)
			{
				if (j != 0) solutionVector.Clear();
				Vector rhsVector = otherMatrix.GetColumn(j);
				factorization.SolveLinearSystem(rhsVector, solutionVector);
				solutionVectors.SetSubcolumn(j, solutionVector);
			}

			watch.Stop();
			Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
			Logger.IncrementAnalysisStep();
			return solutionVectors;
		}

		public class Builder
		{
			public Builder() { }

			public IDofOrderer DofOrderer { get; set; }
				= new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());

			public double FactorizationPivotTolerance { get; set; } = 1E-15;

			public CSparseLUSolver BuildSolver(IModel model)
			{
				return new CSparseLUSolver(model, FactorizationPivotTolerance, DofOrderer);
			}
		}
	}
}
