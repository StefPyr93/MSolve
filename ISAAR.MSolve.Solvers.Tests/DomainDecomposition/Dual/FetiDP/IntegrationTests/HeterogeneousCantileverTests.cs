using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    public static class HeterogeneousCantileverTests
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }
        public enum Residual { Approximate, Exact }

        private const double domainLengthX = 3.0, domainLengthY = 1.5;
        private const int singleSubdomainID = 0;

        //TODO: Exact residual calculation is not implemented yet. Therefore the expected iterations cannot be tested 
        //      accurately using the approximate residual yet.
        [Theory]
        // Homogeneous problem
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.Dirichlet, Residual.Approximate, 11, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.Dirichlet, Residual.Approximate, 11, MatrixFormat.SuiteSparse)]
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.DirichletDiagonal, Residual.Approximate, 14, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.DirichletDiagonal, Residual.Approximate, 14, MatrixFormat.SuiteSparse)]
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.Lumped, Residual.Approximate, 18, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.Homogeneous, 1.0, Precond.Lumped, Residual.Approximate, 18, MatrixFormat.SuiteSparse)]

        // Heterogeneous Lumped
        // Stiffness ratio = 1E-2
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-2, Precond.Dirichlet, Residual.Approximate, 12, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-2, Precond.DirichletDiagonal, Residual.Approximate, 16, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-2, Precond.Lumped, Residual.Approximate, 21, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-3
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-3, Precond.Dirichlet, Residual.Approximate, 13, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-3, Precond.DirichletDiagonal, Residual.Approximate, 20, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-3, Precond.Lumped, Residual.Approximate, 22, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-4
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-4, Precond.Dirichlet, Residual.Approximate, 14, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-4, Precond.DirichletDiagonal, Residual.Approximate, 22, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-4, Precond.Lumped, Residual.Approximate, 26, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-5
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-5, Precond.Dirichlet, Residual.Approximate, 14, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-5, Precond.DirichletDiagonal, Residual.Approximate, 23, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-5, Precond.Lumped, Residual.Approximate, 30, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-6
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-6, Precond.Dirichlet, Residual.Approximate, 15, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-6, Precond.DirichletDiagonal, Residual.Approximate, 27, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousLumped, 1E-6, Precond.Lumped, Residual.Approximate, 33, MatrixFormat.Skyline)]

        // Heterogeneous Condensed
        // Stiffness ratio = 1E-2
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-2, Precond.Dirichlet, Residual.Approximate, 12, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-2, Precond.DirichletDiagonal, Residual.Approximate, 16, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-2, Precond.Lumped, Residual.Approximate, 21, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-3
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-3, Precond.Dirichlet, Residual.Approximate, 13, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-3, Precond.DirichletDiagonal, Residual.Approximate, 20, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-3, Precond.Lumped, Residual.Approximate, 22, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-4
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-4, Precond.Dirichlet, Residual.Approximate, 14, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-4, Precond.DirichletDiagonal, Residual.Approximate, 22, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-4, Precond.Lumped, Residual.Approximate, 26, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-5
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-5, Precond.Dirichlet, Residual.Approximate, 14, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-5, Precond.DirichletDiagonal, Residual.Approximate, 23, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-5, Precond.Lumped, Residual.Approximate, 30, MatrixFormat.Skyline)]
        // Stiffness ratio = 1E-6
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-6, Precond.Dirichlet, Residual.Approximate, 15, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-6, Precond.DirichletDiagonal, Residual.Approximate, 27, MatrixFormat.Skyline)]
        [InlineData(StiffnessDistributionType.HeterogeneousCondensed, 1E-6, Precond.Lumped, Residual.Approximate, 33, MatrixFormat.Skyline)]


        public static void Run(StiffnessDistributionType stiffnessDistributionType, double stiffnessRatio, Precond precond, 
            Residual convergence, int iterExpected, MatrixFormat format)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains(stiffnessRatio);
            (IVectorView ddDisplacements, ISolverLogger logger) = SolveModelWithSubdomains(
                stiffnessDistributionType, stiffnessRatio, precond, convergence, pcgConvergenceTol, format);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            int analysisStep = 0;
            Assert.Equal(140, logger.GetNumDofs(analysisStep, "Lagrange multipliers"));
            Assert.Equal(20, logger.GetNumDofs(analysisStep, "Corner dofs"));

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 6);

            // Allow some tolerance for the iterations:
            int maxIterationsForApproximateResidual = (int)Math.Ceiling(1.0 * iterExpected);
            int pcgIterations = logger.GetNumIterationsOfIterativeAlgorithm(analysisStep);
            Assert.InRange(pcgIterations, 1, maxIterationsForApproximateResidual); // the upper bound is inclusive!
        }

        internal static Model CreateModel(double stiffnessRatio)
        {
            // Subdomains:
            // /|
            // /||-------|-------|-------|-------|  
            // /||  (4)  |  (5)  |  (6)  |  (7)  |
            // /||   E1  |   E0  |   E0  |   E0  |
            // /||-------|-------|-------|-------|  
            // /||  (0)  |  (1)  |  (2)  |  (3)  |
            // /||   E1  |   E0  |   E0  |   E0  |
            // /||-------|-------|-------|-------|
            // /|

            double E0 = 2.1E7;
            double E1 = stiffnessRatio * E0;

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = domainLengthX;
            builder.DomainLengthY = domainLengthY;
            builder.NumSubdomainsX = 4;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 20;
            builder.NumTotalElementsY = 20;
            builder.YoungModuliOfSubdomains = new double[,] { { E1, E0, E0, E0 }, { E1, E0, E0, E0 } };
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 100.0);

            return builder.BuildModel();
        }

        internal static IVectorView SolveModelWithoutSubdomains(double stiffnessRatio)
        {
            Model model = CreateSingleSubdomainModel(stiffnessRatio);

            // Solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[singleSubdomainID].Solution;
        }

        private static Model CreateSingleSubdomainModel(double stiffnessRatio)
        {
            // Replace the existing subdomains with a single one 
            Model model = CreateModel(stiffnessRatio);
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.ElementsDictionary.Values) subdomain.Elements.Add(element.ID, element);
            return model;
        }

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(
            StiffnessDistributionType stiffnessDistributionType, double stiffnessRatio,
            Precond precond, Residual residualConvergence, double pcgConvergenceTolerance, MatrixFormat format)
        {
            // Model
            Model multiSubdomainModel = CreateModel(stiffnessRatio);
            multiSubdomainModel.ConnectDataStructures();

            // Corner nodes
            double meshTol = 1E-6;
            var cornerNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in multiSubdomainModel.SubdomainsDictionary.Values)
            {
                subdomain.DefineNodesFromElements(); //TODO: This will also be called by the analyzer.
                INode[] corners = CornerNodeUtilities.FindCornersOfRectangle2D(subdomain);
                var cornerNodes = new HashSet<INode>();
                foreach (INode node in corners)
                {
                    if (node.Constraints.Count > 0) continue;
                    if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y) <= meshTol)) continue;
                    if ((Math.Abs(node.X - domainLengthX) <= meshTol) && (Math.Abs(node.Y - domainLengthY) <= meshTol)) continue;
                    cornerNodes.Add(node);
                }
                cornerNodesOfEachSubdomain[subdomain] = cornerNodes;
            }
            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodesOfEachSubdomain);


            // Solver
            IFetiDPMatrixManagerFactory fetiMatrices = MatrixFormatSelection.DefineMatrixManagerFactory(format);
            var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
            solverBuilder.StiffnessDistribution = stiffnessDistributionType;


            // Preconditioner
            if (precond == Precond.Lumped) solverBuilder.Preconditioning = new LumpedPreconditioning();
            else if (precond == Precond.Dirichlet) solverBuilder.Preconditioning = new DirichletPreconditioning();
            else solverBuilder.Preconditioning = new DiagonalDirichletPreconditioning();

            //TODO: This needs to be implemented for FETI-DP
            //// PCG may need to use the exact residual for the comparison with the expected values
            //bool residualIsExact = residualConvergence == Residual.Exact;
            //ExactPcpgConvergence.Factory exactResidualConvergence = null;
            //if (residualIsExact)
            //{
            //    exactResidualConvergence = new ExactPcpgConvergence.Factory(
            //        CreateSingleSubdomainModel(stiffnessRatio), solverBuilder.DofOrderer,
            //            (model, solver) => new ProblemStructural(model, solver));
            //}
            //if (residualIsExact) exactResidualConvergence.InterfaceProblemSolver = interfaceProblemSolver;

            // Specify PCG settings
            solverBuilder.PcgSettings = new PcgSettings() { ConvergenceTolerance = pcgConvergenceTolerance };

            FetiDPSolverSerial fetiSolver = solverBuilder.Build(multiSubdomainModel, cornerNodeSelection);
            //if (residualIsExact) exactResidualConvergence.FetiSolver = fetiSolver;

            // Run the analysis
            FetiDPSolverSerialTests.RunAnalysis(multiSubdomainModel, fetiSolver); //check dof separator

            // Gather the global displacements
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();

            return (globalDisplacements, fetiSolver.Logger);
        }
    }
}
