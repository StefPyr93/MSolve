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
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    public static class Cantilever2DFetiDPTests
    {
        public enum Corners { C2424, C0202, C1313, C0213 }
        public enum Crosspoints { Minimum, FullyRedundant }
        public enum Augmented { None, Node21}

        private const int singleSubdomainID = 0;
        
        [Theory]
        // Homogeneous problem
        [InlineData(Corners.C2424, Crosspoints.Minimum, Augmented.None)]
        [InlineData(Corners.C1313, Crosspoints.Minimum, Augmented.None)]
        [InlineData(Corners.C1313, Crosspoints.FullyRedundant, Augmented.None)]
        [InlineData(Corners.C1313, Crosspoints.Minimum, Augmented.Node21)]
        [InlineData(Corners.C1313, Crosspoints.FullyRedundant, Augmented.Node21)]
        //[InlineData(Corners.C0213, Crosspoints.Minimum, Augmented.None)]
        //[InlineData(Corners.C0213, Crosspoints.FullyRedundant, Augmented.None)]
        //[InlineData(Corners.C0202, Crosspoints.Minimum, Augmented.None)]
        //[InlineData(Corners.C0202, Crosspoints.FullyRedundant, Augmented.None)]
        public static void Run(Corners corners, Crosspoints crosspoints, Augmented augmented)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains();
            (IVectorView ddDisplacements, ISolverLogger logger) =
                SolveModelWithSubdomains(corners, crosspoints, augmented, pcgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 5);
        }

        internal static Model CreateModel()
        {
            // Subdomains:
            // /|
            // /||-------|-------| |-------|-------|  
            // /||       |       | |       |       |
            // /||       |       | |       |       |
            // /||-------|-------| |-------|-------|
            //                     
            // /||-------|-------| |-------|-------|  
            // /||       |       | |       |       |
            // /||       |       | |       |       |
            // /||-------|-------| |-------|-------|
            // /|

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 4.0;
            builder.DomainLengthY = 2.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 4;
            builder.NumTotalElementsY = 2;
            builder.YoungModulus = 2.1E7;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, -600.0);

            return builder.BuildModel();
        }

        internal static IVectorView SolveModelWithoutSubdomains()
        {
            Model model = CreateSingleSubdomainModel();

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

        private static Model CreateSingleSubdomainModel()
        {
            // Replace the existing subdomains with a single one 
            Model model = CreateModel();
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.ElementsDictionary.Values) subdomain.Elements.Add(element.ID, element);
            return model;
        }

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(Corners corners,
            Crosspoints crosspoints, Augmented augmented, double pcgConvergenceTolerance)
        {
            // Model
            Model model = CreateModel();
            model.ConnectDataStructures();

            // Corner, midside nodes
            ICornerNodeSelection cornerNodeSelection = DefineCornerNodes(model, corners);
            IMidsideNodesSelection midsideNodes = DefineMidsideNodes(model, augmented);

            // Crosspoint strategy
            ICrosspointStrategy crosspointStrategy = null;
            if (crosspoints == Crosspoints.FullyRedundant) crosspointStrategy = new FullyRedundantConstraints();
            else if (crosspoints == Crosspoints.Minimum) crosspointStrategy = new MinimumConstraints();
            else throw new ArgumentException();

            // Specify PCG settings
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = pcgConvergenceTolerance,
                MaxIterationsProvider = new FixedMaxIterationsProvider(100)
            };

            // Solver
            if (augmented == Augmented.None)
            {
                var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
                var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
                solverBuilder.Preconditioning = new DirichletPreconditioning();
                solverBuilder.CrosspointStrategy = crosspointStrategy;
                solverBuilder.PcgSettings = pcgSettings;
                FetiDPSolverSerial fetiSolver = solverBuilder.Build(model, cornerNodeSelection);

                // Run the analysis
                RunAnalysis(model, fetiSolver); //check dof separator

                // Gather the global displacements
                Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();
                return (globalDisplacements, fetiSolver.Logger);
            }
            else if (augmented == Augmented.Node21)
            {

                var fetiMatrices = new FetiDP3dMatrixManagerFactoryDense();
                var solverBuilder = new FetiDP3dSolverSerial.Builder(fetiMatrices);
                solverBuilder.Preconditioning = new DirichletPreconditioning();
                solverBuilder.CrosspointStrategy = crosspointStrategy;
                solverBuilder.PcgSettings = pcgSettings;
                FetiDP3dSolverSerial fetiSolver = solverBuilder.Build(model, cornerNodeSelection, midsideNodes);

                // Run the analysis
                RunAnalysis(model, fetiSolver); //check dof separator

                // Gather the global displacements
                Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();
                return (globalDisplacements, fetiSolver.Logger);
            }
            else throw new ArgumentException();
        }

        private static ICornerNodeSelection DefineCornerNodes(Model model, Corners corners)
        {
            // Important nodes
            double meshTol = 1E-6;
            Node n20 = FindNode(2, 0, model, meshTol);
            Node n21 = FindNode(2, 1, model, meshTol);
            Node n22 = FindNode(2, 2, model, meshTol);
            IEnumerable<Node> nodes2 = FindNodesWithX(2, model, meshTol);

            Node n40 = FindNode(4, 0, model, meshTol);
            Node n41 = FindNode(4, 1, model, meshTol);
            Node n42 = FindNode(4, 2, model, meshTol);
            IEnumerable<Node> nodes4 = FindNodesWithX(4, model, meshTol);

            // Corner nodes
            var cornerNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                cornerNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            if (corners == Corners.C2424)
            {
                foreach (Node node in nodes4.Union(nodes2))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C0202)
            {
                foreach (Node node in nodes4)
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C1313)
            {
                var extraCorners = new Node[] { n20, n22 };
                foreach (Node node in nodes4.Union(extraCorners))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C0213)
            {
                var extraCorners = new Node[] { n22 };
                foreach (Node node in nodes4.Union(extraCorners))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else throw new ArgumentException();
            return new UsedDefinedCornerNodes(cornerNodesOfEachSubdomain);
        }

        private static IMidsideNodesSelection DefineMidsideNodes(Model model, Augmented augmented)
        {
            // Important nodes
            double meshTol = 1E-6;
            Node n20 = FindNode(2, 0, model, meshTol);
            Node n21 = FindNode(2, 1, model, meshTol);
            Node n22 = FindNode(2, 2, model, meshTol);
            IEnumerable<Node> nodes2 = FindNodesWithX(2, model, meshTol);

            Node n40 = FindNode(4, 0, model, meshTol);
            Node n41 = FindNode(4, 1, model, meshTol);
            Node n42 = FindNode(4, 2, model, meshTol);
            IEnumerable<Node> nodes4 = FindNodesWithX(4, model, meshTol);

            // Midside nodes
            var midsideNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                midsideNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            if (augmented == Augmented.None) return null;
            else if (augmented == Augmented.Node21)
            {
                foreach (Subdomain subdomain in n21.SubdomainsDictionary.Values)
                {
                    midsideNodesOfEachSubdomain[subdomain].Add(n21);
                }
            }
            else throw new ArgumentException();
            return new UserDefinedMidsideNodes(midsideNodesOfEachSubdomain, 
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY});
        }

        private static Node FindNode(double x, double y, Model model, double tol)
        {
            return model.NodesDictionary.Values.Where(
                n => (Math.Abs(n.X - x) <= tol) && (Math.Abs(n.Y - y) <= tol)).First();
        }

        private static IEnumerable<Node> FindNodesWithX(double x, Model model, double tol)
        {
            return model.NodesDictionary.Values.Where(n => (Math.Abs(n.X - x) <= tol));
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
