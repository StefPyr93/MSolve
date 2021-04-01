using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    public static class Rve2DFetiDPTests
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }
        public enum Crosspoints { Minimum, FullyRedundant }

        private const double domainLengthX = 3.0, domainLengthY = 1.5;
        private const int singleSubdomainID = 0;
        
        [Theory]
        // Homogeneous problem
        [InlineData(1.0, Precond.Dirichlet, Crosspoints.FullyRedundant)]
        [InlineData(1.0, Precond.DirichletDiagonal, Crosspoints.FullyRedundant)]
        [InlineData(1.0, Precond.Lumped, Crosspoints.FullyRedundant)]
        [InlineData(1.0, Precond.Dirichlet, Crosspoints.Minimum)]
        [InlineData(1.0, Precond.DirichletDiagonal, Crosspoints.Minimum)]
        [InlineData(1.0, Precond.Lumped, Crosspoints.Minimum)]
        public static void Run(double stiffnessRatio, Precond precond, Crosspoints crosspoints)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains(stiffnessRatio);
            (IVectorView ddDisplacements, ISolverLogger logger) =
                SolveModelWithSubdomains(stiffnessRatio, precond, crosspoints, pcgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 5);
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
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.UpperSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.UpperSide, StructuralDof.TranslationY, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerSide, StructuralDof.TranslationY, 0.0);

            Model model = builder.BuildModel();
            double tol = 1E-6;
            IEnumerable<Node> internalNodes = model.NodesDictionary.Values.Where(node =>
                (Math.Abs(node.X - 0.5 * domainLengthX) <= tol) && (Math.Abs(node.Y - 0.5 * domainLengthY) <= tol));
            Debug.Assert(internalNodes.Count() >= 1);
            model.Loads.Add(new Load() { Node = internalNodes.First(), DOF = StructuralDof.TranslationY, Amount = 1000 });

            return model;
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

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(double stiffnessRatio,
            Precond precond, Crosspoints crosspoints, double pcgConvergenceTolerance)
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
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
            var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
            if (stiffnessRatio == 1.0) solverBuilder.StiffnessDistribution = StiffnessDistributionType.Homogeneous;
            else solverBuilder.StiffnessDistribution = StiffnessDistributionType.HeterogeneousLumped;

            // Preconditioner
            if (precond == Precond.Lumped) solverBuilder.Preconditioning = new LumpedPreconditioning();
            else if (precond == Precond.Dirichlet) solverBuilder.Preconditioning = new DirichletPreconditioning();
            else solverBuilder.Preconditioning = new DiagonalDirichletPreconditioning();

            // Crosspoint strategy
            if (crosspoints == Crosspoints.FullyRedundant) solverBuilder.CrosspointStrategy = new FullyRedundantConstraints();
            else if (crosspoints == Crosspoints.Minimum) solverBuilder.CrosspointStrategy = new MinimumConstraints();
            else throw new ArgumentException();

            // Specify PCG settings
            solverBuilder.PcgSettings = new PcgSettings() { ConvergenceTolerance = pcgConvergenceTolerance };

            FetiDPSolverSerial fetiSolver = solverBuilder.Build(multiSubdomainModel, cornerNodeSelection);
            //if (residualIsExact) exactResidualConvergence.FetiSolver = fetiSolver;

            // Run the analysis
            RunAnalysis(multiSubdomainModel, fetiSolver); //check dof separator

            // Gather the global displacements
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();

            return (globalDisplacements, fetiSolver.Logger);
        }

        private static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis
            solver.OrderDofs(false);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
            solver.BuildGlobalMatrix(new ElementStructuralStiffnessProvider());
            model.ApplyLoads();
            LoadingUtilities.ApplyNodalLoads(model, solver);
            solver.Solve();
        }
    }
}
