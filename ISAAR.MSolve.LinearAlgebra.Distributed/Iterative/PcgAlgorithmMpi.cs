using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using MPI;

//TODO: Remove/reduce duplication between this and the serial version.
//TODO: perhaps all quantities should be stored as mutable fields, exposed as readonly properties and the various strategies 
//      should read them from a reference of CG/PCG/PCPG, instead of having them injected.
//TODO: In regular CG, there is a check to perevent premature convergence, by correcting the residual. Can this be done for PCG 
//      as well? Would the preconditioned residual be updated as well?
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Iterative
{
    /// <summary>
    /// Implements the untransformed Preconditioned Conjugate Gradient algorithm for solving linear systems with symmetric 
    /// positive definite matrices. This implementation is based on the algorithm presented in section B3 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgAlgorithmMpi : PcgAlgorithmBase
    {
        private const string name = "Preconditioned Conjugate Gradient";
        private readonly IPcgBetaParameterCalculation betaCalculation;
        private readonly Intracommunicator comm; //TODO: Use a dedicated class to handle this and master.
        private readonly int master;

        private PcgAlgorithmMpi(Intracommunicator comm, int masterProcess, double residualTolerance,
            IMaxIterationsProvider maxIterationsProvider, IPcgResidualConvergence pcgConvergence,
            IPcgResidualUpdater residualUpdater, IPcgBetaParameterCalculation betaCalculation) : 
            base(residualTolerance, maxIterationsProvider, pcgConvergence, residualUpdater)
        {
            this.comm = comm;
            this.master = masterProcess;
            this.betaCalculation = betaCalculation;
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">
        /// Represents the matrix A of the linear system A * x = b, which must be symmetric positive definite.
        /// </param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="preconditioner">
        /// A preconditioner matrix that is also symmetric positive definite and has the same dimensions as A.
        /// </param>
        /// <param name="solution">
        /// The vector from which to start refining the solution vector x. Constraints:
        /// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="initialGuessIsZero">
        /// If <paramref name="solution"/> is 0, then set <paramref name="initialGuessIsZero"/> to true to avoid performing the
        /// operation b-A*0 before starting.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        public override IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
            IVector solution, bool initialGuessIsZero, Func<IVector> zeroVectorInitializer)
        {
            if (comm.Rank == master)
            {
                //TODO: find a better way to handle optimizations for the case x0=0, than using an initialGuessIsZero flag
                Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
                Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);
            }

            this.Matrix = matrix;
            this.Preconditioner = preconditioner;
            this.Rhs = rhs;
            this.solution = solution;

            // r = b - A * x
            if (comm.Rank == master)
            {
                if (initialGuessIsZero) residual = rhs.Copy();
                else residual = ExactResidual.Calculate(matrix, rhs, solution);
            }

            return SolveInternal(maxIterationsProvider.GetMaxIterations(matrix.NumColumns), zeroVectorInitializer);
        }

        protected override IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer)
        {
            if (comm.Rank == master)
            {
                // In contrast to the source algorithm, we initialize s here. At each iteration it will be overwritten, 
                // thus avoiding allocating & deallocating a new vector.
                precondResidual = zeroVectorInitializer();

                // d = inv(M) * r
                direction = zeroVectorInitializer();
            }
            Preconditioner.SolveLinearSystem(residual, direction);

            if (comm.Rank == master)
            {
                // δnew = δ0 = r * d
                resDotPrecondRes = residual.DotProduct(direction);

                // The convergence and beta strategies must be initialized immediately after the first r and r*inv(M)*r are computed.
                convergence.Initialize(this);
                betaCalculation.Initialize(this);

                // Allocate memory for other vectors, which will be reused during each iteration
                matrixTimesDirection = zeroVectorInitializer();
            }

            // This is also used as output
            double residualNormRatio = double.NaN;

            for (iteration = 0; iteration < maxIterations; ++iteration)
            {
                // q = A * d
                Matrix.Multiply(direction, matrixTimesDirection);

                if (comm.Rank == master)
                {
                    // α = δnew / (d * q)
                    stepSize = resDotPrecondRes / direction.DotProduct(matrixTimesDirection);

                    // x = x + α * d
                    solution.AxpyIntoThis(direction, stepSize);

                    // Normally the residual vector is updated as: r = r - α * q. However corrections might need to be applied.
                    residualUpdater.UpdateResidual(this, residual);
                }

                // s = inv(M) * r
                Preconditioner.SolveLinearSystem(residual, precondResidual);

                if (comm.Rank == master)
                {
                    // δold = δnew
                    resDotPrecondResOld = resDotPrecondRes;

                    // δnew = r * s 
                    resDotPrecondRes = residual.DotProduct(precondResidual);

                    // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                    residualNormRatio = convergence.EstimateResidualNormRatio(this);
                }

                comm.Broadcast<double>(ref residualNormRatio, master);
                Debug.WriteLine($"Process {comm.Rank}: PCG Iteration = {iteration}: residual norm ratio = {residualNormRatio}");
                if (residualNormRatio <= residualTolerance)
                {
                    return new IterativeStatistics
                    {
                        AlgorithmName = name,
                        HasConverged = true,
                        NumIterationsRequired = iteration + 1,
                        ResidualNormRatioEstimation = residualNormRatio
                    };
                }

                if (comm.Rank == master)
                {
                    // The default Fletcher-Reeves formula is: β = δnew / δold = (sNew * rNew) / (sOld * rOld)
                    // However we could use a different one, e.g. for variable preconditioning Polak-Ribiere is usually better.
                    paramBeta = betaCalculation.CalculateBeta(this);

                    // d = s + β * d
                    //TODO: benchmark the two options to find out which is faster
                    //direction = preconditionedResidual.Axpy(direction, beta); //This allocates a new vector d, copies r and GCs the existing d.
                    direction.LinearCombinationIntoThis(paramBeta, precondResidual, 1.0); //This performs additions instead of copying and needless multiplications.
                }
            }

            // We reached the max iterations before PCG converged
            return new IterativeStatistics
            {
                AlgorithmName = name,
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                ResidualNormRatioEstimation = residualNormRatio
            };
        }

        /// <summary>
        /// Constructs <see cref="PcgAlgorithm"/> instances, allows the user to specify some or all of the required parameters 
        /// and provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder : PcgBuilderBase
        {
            /// <summary>
            /// Specifies how to calculate the beta parameter of PCG, which is used to update the direction vector. 
            /// </summary>
            public IPcgBetaParameterCalculation BetaCalculation { get; set; } = new FletcherReevesBeta();

            /// <summary>
            /// Creates a new instance of <see cref="PcgAlgorithmMpi"/>.
            /// </summary>
            public PcgAlgorithmMpi Build(Intracommunicator comm, int masterProcess)
            {
                return new PcgAlgorithmMpi(comm, masterProcess, ResidualTolerance, MaxIterationsProvider, Convergence,
                    ResidualUpdater, BetaCalculation);
            }
        }
    }
}
