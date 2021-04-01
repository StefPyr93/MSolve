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
using ISAAR.MSolve.Logging.DomainDecomposition;
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
    public static class Cantilever3DFetiDPTests
    {
        private const string plotPath = @"C:\Users\Serafeim\Desktop\FETI-DP\Plots";
        private const double meshTol = 1E-6;

        private const int singleSubdomainID = 0;

        public enum Crosspoints { Minimum, FullyRedundant }
        public enum Corners 
        {
            Case0, Case1, Case2, Case3,
            C0426, C0404, C1526
        }
        public enum Augmentation { None, Case0, Case1 }
        public enum Mesh { e4x2x1, e4x4x2 }

        //[Theory]
        //[InlineData(Mesh.e4x2x1, Corners.Case0, Crosspoints.Minimum, Augmentation.None)]
        //[InlineData(Mesh.e4x2x1, Corners.Case1, Crosspoints.Minimum, Augmentation.None)]
        //[InlineData(Mesh.e4x2x1, Corners.Case1, Crosspoints.FullyRedundant, Augmentation.None)]
        //[InlineData(Mesh.e4x2x1, Corners.Case1, Crosspoints.Minimum, Augmentation.Case0)]
        //[InlineData(Mesh.e4x2x1, Corners.Case1, Crosspoints.FullyRedundant, Augmentation.Case0)]

        //[InlineData(Mesh.e4x4x2, Corners.Case2, Crosspoints.Minimum, Augmentation.None)]
        //[InlineData(Mesh.e4x4x2, Corners.Case2, Crosspoints.FullyRedundant, Augmentation.None)]
        //[InlineData(Mesh.e4x4x2, Corners.Case2, Crosspoints.Minimum, Augmentation.Case1)]
        //[InlineData(Mesh.e4x4x2, Corners.Case2, Crosspoints.FullyRedundant, Augmentation.Case1)]
        //[InlineData(Mesh.e4x4x2, Corners.Case3, Crosspoints.Minimum, Augmentation.None)]
        //[InlineData(Mesh.e4x4x2, Corners.Case3, Crosspoints.FullyRedundant, Augmentation.None)]
        //[InlineData(Mesh.e4x4x2, Corners.Case3, Crosspoints.Minimum, Augmentation.Case1)]
        //[InlineData(Mesh.e4x4x2, Corners.Case3, Crosspoints.FullyRedundant, Augmentation.Case1)]

        //[InlineData(Corners.C1526, Crosspoints.Minimum, Augmented.None)]
        //[InlineData(Corners.C1526, Crosspoints.FullyRedundant, Augmented.None)]
        // Each subdomain must have at least 1 corner node. Otherwise some matrices degenerate
        //[InlineData(Corners.C0426, Crosspoints.Minimum, Augmented.None)]          
        //[InlineData(Corners.C0426, Crosspoints.FullyRedundant, Augmented.None)]
        //[InlineData(Corners.C0404, Crosspoints.Minimum)]
        //[InlineData(Corners.C0404, Crosspoints.FullyRedundant)]
        public static void Run(Mesh mesh, Corners corners, Crosspoints crosspoints, Augmentation augmentation)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains(CreateModel(mesh));
            (IVectorView ddDisplacements, ISolverLogger logger) =
                SolveModelWithSubdomains(CreateModel(mesh), corners, crosspoints, augmentation, pcgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 6);
            //Assert.True(directDisplacements.Equals(ddDisplacements, 1E-5));
        }

        internal static Model CreateModel(Mesh mesh)
        {
            int numElementsX = -1, numElementsY = -1, numElementsZ = -1;
            if (mesh == Mesh.e4x2x1)
            {
                numElementsX = 4; numElementsY = 2; numElementsZ = 1;
            }
            else if (mesh == Mesh.e4x4x2)
            {
                numElementsX = 4; numElementsY = 4; numElementsZ = 2;
            }

            double E0 = 2.1E7;
            //double E1 = stiffnessRatio * E0;

            var builder = new Uniform3DModelBuilder();
            builder.MaxX = 4.0;
            builder.MaxY = 2.0;
            builder.MaxZ = 1.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumSubdomainsZ = 1;
            builder.NumTotalElementsX = numElementsX;
            builder.NumTotalElementsY = numElementsY;
            builder.NumTotalElementsZ = numElementsZ;
            builder.YoungModulus = E0;
            //builder.YoungModuliOfSubdomains = new double[,] { { E1, E0, E0, E0 }, { E1, E0, E0, E0 } };

            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationY, 0.0);
            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationZ, 0.0);
            builder.DistributeLoadAtNodes(Uniform3DModelBuilder.BoundaryRegion.MaxX, StructuralDof.TranslationY, -600.0);

            return builder.BuildModel();
        }

        internal static IVectorView SolveModelWithoutSubdomains(Model model)
        {
            // Replace the existing subdomains with a single one 
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.ElementsDictionary.Values) subdomain.Elements.Add(element.ID, element);

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

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(Model model,
            Corners corners, Crosspoints crosspoints, Augmentation augmentation, double pcgConvergenceTolerance)
        {
            // Model
            model.ConnectDataStructures();

            // Corner, midside nodes
            ICornerNodeSelection cornerNodeSelection = DefineCornerNodes(model, corners);
            IMidsideNodesSelection midsideNodes = DefineMidsideNodes(model, augmentation);

            // Crosspoint strategy
            ICrosspointStrategy crosspointStrategy = null;
            if (crosspoints == Crosspoints.FullyRedundant) crosspointStrategy = new FullyRedundantConstraints();
            else if (crosspoints == Crosspoints.Minimum) crosspointStrategy = new MinimumConstraints();
            else throw new ArgumentException();

            //// Plot for debugging
            //var logger = new DomainDecompositionLoggerFetiDP(cornerNodeSelection, plotPath, true);
            //logger.PlotSubdomains(model);

            // Specify PCG settings
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = pcgConvergenceTolerance,
                MaxIterationsProvider = new FixedMaxIterationsProvider(100)
            };

            // Solver
            if (augmentation == Augmentation.None)
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
            else
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
        }

        private static ICornerNodeSelection DefineCornerNodes(Model model, Corners corners)
        {
            // Select corner nodes
            var cornerNodes = new HashSet<Node>();
            if (corners == Corners.Case0)
            {
                cornerNodes.UnionWith(FindNodesWithX(4, model));
                cornerNodes.UnionWith(FindNodesWithX(2, model));
            }
            else if (corners == Corners.Case1)
            {
                cornerNodes.UnionWith(FindNodesWithX(4, model));
                cornerNodes.Add(FindNode(2, 0, 0, model));
                cornerNodes.Add(FindNode(2, 2, 0, model));
                cornerNodes.Add(FindNode(2, 0, 1, model));
                cornerNodes.Add(FindNode(2, 2, 1, model));
            }
            else if (corners == Corners.Case2)
            {
                cornerNodes.Add(FindNode(2, 1, 0, model));
                cornerNodes.Add(FindNode(2, 1, 1, model));
                cornerNodes.Add(FindNode(4, 1, 1, model));
            }
            else if (corners == Corners.Case3)
            {
                cornerNodes.Add(FindNode(2, 1, 0, model));
                cornerNodes.Add(FindNode(2, 1, 1, model));
                cornerNodes.Add(FindNode(4, 1, 1, model));
                cornerNodes.Add(FindNode(2, 0, 1, model));
                cornerNodes.Add(FindNode(2, 2, 1, model));
            }
            else if (corners == Corners.C0404)
            {
                cornerNodes.UnionWith(FindNodesWithX(4, model));
            }
            else if (corners == Corners.C0426)
            {
                cornerNodes.UnionWith(FindNodesWithX(4, model));
                cornerNodes.Add(FindNode(2, 2, 0, model));
                cornerNodes.Add(FindNode(2, 2, 1, model));
            }
            else if (corners == Corners.C1526)
            {
                cornerNodes.UnionWith(FindNodesWithX(4, model));
                cornerNodes.Add(FindNode(2, 2, 0, model));
                cornerNodes.Add(FindNode(2, 0, 1, model));
                cornerNodes.Add(FindNode(2, 2, 1, model));
            }
            else throw new ArgumentException();

            var cornerNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                cornerNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            foreach (Node node in cornerNodes)
            {
                foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                {
                    cornerNodesOfEachSubdomain[subdomain].Add(node);
                }
            }
            return new UsedDefinedCornerNodes(cornerNodesOfEachSubdomain);
        }

        private static IMidsideNodesSelection DefineMidsideNodes(Model model, Augmentation augmentation)
        {
            // Select midside nodes
            var midsideNodes = new HashSet<Node>();
            if (augmentation == Augmentation.None) return null;
            else if (augmentation == Augmentation.Case0)
            {
                midsideNodes.Add(FindNode(2, 1, 0, model));
                midsideNodes.Add(FindNode(2, 1, 1, model));
            }
            else if (augmentation == Augmentation.Case1)
            {
                midsideNodes.Add(FindNode(2, 1, 0.5, model));
            }
            else throw new ArgumentException();

            // Midside nodes
            var midsideNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                midsideNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            foreach (Node node in midsideNodes)
            {
                foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                {
                    midsideNodesOfEachSubdomain[subdomain].Add(node);
                }
            }
            return new UserDefinedMidsideNodes(midsideNodesOfEachSubdomain,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
        }

        private static Node FindNode(double x, double y, double z, Model model)
        {
            return model.NodesDictionary.Values.Where(n => 
                (Math.Abs(n.X - x) <= meshTol) && (Math.Abs(n.Y - y) <= meshTol) && (Math.Abs(n.Z - z) <= meshTol)).First();
        }

        private static IEnumerable<Node> FindNodesWithX(double x, Model model)
        {
            return model.NodesDictionary.Values.Where(n => (Math.Abs(n.X - x) <= meshTol));
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
