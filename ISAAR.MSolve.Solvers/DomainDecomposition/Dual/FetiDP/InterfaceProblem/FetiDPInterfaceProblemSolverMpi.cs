using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;
using MPI;
using ISAAR.MSolve.LinearAlgebra.Distributed.Iterative;
using ISAAR.MSolve.LinearAlgebra.Distributed.Vectors;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;

//TODO: Reduce the duplication between MPI and serial implementations. Most is FETI-DP specific code.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// The interface problem is solved using PCG. The matrix of the coarse problem KccStar, namely the static condensation of 
    /// the remainder dofs onto the corner dofs is performed explicitly.
    /// </summary>
    public class FetiDPInterfaceProblemSolverMpi : IFetiDPInterfaceProblemSolver
    {
        private readonly IModel model;
        private readonly PcgSettings pcgSettings;
        private readonly ProcessDistribution procs;

        public FetiDPInterfaceProblemSolverMpi(ProcessDistribution processDistribution, IModel model, PcgSettings pcgSettings)
        {
            this.procs = processDistribution;
            this.model = model;
            this.pcgSettings = pcgSettings;
        }

        public Vector PreviousLambda { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public bool UsePreviousLambda { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public ReorthogonalizedPcg Pcg => throw new NotImplementedException();

        public bool UseStagnationCriterion { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public Vector SolveInterfaceProblem(IFetiDPMatrixManager matrixManager,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, double globalForcesNorm, ISolverLogger logger)
        {
            int systemOrder = flexibility.NumGlobalLagrangeMultipliers;

            // Prepare PCG matrix, preconditioner, rhs and solution
            var pcgMatrix = new FetiDPInterfaceProblemMatrixMpi(procs, matrixManager, flexibility);
            var pcgPreconditioner = new FetiDPInterfaceProblemPreconditioner(preconditioner);
            Vector pcgRhs = null;
            Vector lagranges = null;
            Vector globalDr = CalcGlobalDr(matrixManager, lagrangesEnumerator);
            pcgRhs = CalcInterfaceProblemRhs(matrixManager, flexibility, globalDr);
            if (procs.IsMasterProcess) lagranges = Vector.CreateZero(systemOrder);

            // Solve the interface problem using PCG algorithm
            var pcgBuilder = new PcgAlgorithmMpi.Builder();
            pcgBuilder.MaxIterationsProvider = pcgSettings.MaxIterationsProvider;
            pcgBuilder.ResidualTolerance = pcgSettings.ConvergenceTolerance;
            pcgBuilder.Convergence = pcgSettings.ConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            PcgAlgorithmMpi pcg = pcgBuilder.Build(procs.Communicator, procs.MasterProcess); //TODO: perhaps use the pcg from the previous analysis if it has reorthogonalization.

            IterativeStatistics stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagranges, true,
                () => Vector.CreateZero(systemOrder));

            // Log statistics about PCG execution
            FetiDPInterfaceProblemUtilities.CheckConvergence(stats);
            logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);

            return lagranges;
        }

        private Vector CalcGlobalDr(IFetiDPMatrixManager matrixManager, ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            var transferrer = new VectorTransferrer(procs);
            var subdomainContributions = new List<Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomainContributions.Add(
                    FetiDPInterfaceProblemUtilities.CalcSubdomainDr(subdomain, matrixManager, lagrangesEnumerator));
            }
            return transferrer.SumVectors(subdomainContributions);
        }

        private Vector CalcInterfaceProblemRhs(IFetiDPMatrixManager matrixManager, IFetiDPFlexibilityMatrix flexibility,
            Vector globalDr)
        {
            // rhs = dr - FIrc * inv(KccStar) * fcStar
            Vector temp = null;
            if (procs.IsMasterProcess) temp = matrixManager.MultiplyInverseCoarseProblemMatrix(matrixManager.CoarseProblemRhs);
            temp = flexibility.MultiplyFIrc(temp);
            if (procs.IsMasterProcess) return globalDr - temp;
            else return null;
        }
    }
}
