using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using MPI;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    /// <summary>
    /// Tests from Papagiannakis bachelor thesis (NTUA 2011), p. 134 - 147
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class PapagiannakisFetiDPTests2DMpi
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }
        public enum Residual { Approximate, Exact }

        private const double domainLengthX = 3.0, domainLengthY = 1.5;

        public static void Run(int numProcesses, double stiffnessRatio, Precond precond, Residual convergence, int iterExpected, 
            MatrixFormat format)
        {
            var procs = ProcessDistribution.CreateDistribution(numProcesses, 8);

            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = null;
            if (procs.IsMasterProcess)
            {
                directDisplacements = PapagiannakisFetiDPTests2DSerial.SolveModelWithoutSubdomains(stiffnessRatio);
            }

            (IVectorView ddDisplacements, ISolverLogger logger) =
                SolveModelWithSubdomains(procs, stiffnessRatio, precond, convergence, pcgConvergenceTol, format);

            //Checks
            if (procs.IsMasterProcess)
            {
                int analysisStep = 0;
                Assert.Equal(140, logger.GetNumDofs(analysisStep, "Lagrange multipliers"));
                Assert.Equal(20, logger.GetNumDofs(analysisStep, "Corner dofs"));

                // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
                // different machines to achieve the exact same accuracy.

                double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();
                Assert.Equal(0.0, normalizedError, 6);

                // Allow some tolerance for the iterations:
                int maxIterationsForApproximateResidual = (int)Math.Ceiling(1.0 * iterExpected);
                int pcgIterations = logger.GetNumIterationsOfIterativeAlgorithm(analysisStep);
                Assert.InRange(pcgIterations, 1, maxIterationsForApproximateResidual); // the upper bound is inclusive!
            }
        }

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(ProcessDistribution procs,
            double stiffnessRatio, Precond precond, Residual residualConvergence, double pcgConvergenceTolerance, 
            MatrixFormat format)
        {
            // Model
            //var model = new ModelMpiCentralized(procs, () => PapagiannakisFetiDPTests2DSerial.CreateModel(stiffnessRatio));
            var model = new ModelMpiRedundant(procs, () => PapagiannakisFetiDPTests2DSerial.CreateModel(stiffnessRatio));
            model.ConnectDataStructures();
            model.ScatterSubdomains();

            // Corner nodes
            var cornerNodes = new Dictionary<ISubdomain, HashSet<INode>>();
            if (procs.IsMasterProcess)
            {
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    cornerNodes[subdomain] = DefineCornerNodes(subdomain);
                }
            }
            else
            {
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    ISubdomain subdomain = model.GetSubdomain(s);
                    cornerNodes[subdomain] = DefineCornerNodes(subdomain);
                }
            }
            var cornerNodeSelection = new UsedDefinedCornerNodesMpi(procs, cornerNodes);

            // Solver
            IFetiDPMatrixManagerFactory fetiMatrices = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            var solverBuilder = new FetiDPSolverMpi.Builder(procs, fetiMatrices);
            solverBuilder.ProblemIsHomogeneous = stiffnessRatio == 1.0;
            solverBuilder.PcgSettings = new PcgSettings() { ConvergenceTolerance = pcgConvergenceTolerance };

            // Preconditioner
            if (precond == Precond.Lumped) solverBuilder.Preconditioning = new LumpedPreconditioning();
            else if (precond == Precond.Dirichlet) solverBuilder.Preconditioning = new DirichletPreconditioning();
            else solverBuilder.Preconditioning = new DiagonalDirichletPreconditioning();

            FetiDPSolverMpi fetiSolver = solverBuilder.Build(model, cornerNodeSelection);

            // Run the analysis
            FetiDPSolverMpiTests.RunAnalysis(procs, model, fetiSolver); //check dof separator

            // Gather the global displacements
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();

            return (globalDisplacements, fetiSolver.Logger);
        }

        private static HashSet<INode> DefineCornerNodes(ISubdomain subdomain)
        {
            double meshTol = 1E-6;
            INode[] corners = CornerNodeUtilities.FindCornersOfRectangle2D(subdomain);
            var cornerNodes = new HashSet<INode>();
            foreach (INode node in corners)
            {
                if (node.Constraints.Count > 0) continue;
                if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y) <= meshTol)) continue;
                if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y - domainLengthY) <= meshTol)) continue;
                cornerNodes.Add(node);
            }
            return cornerNodes;
        }
    }
}
