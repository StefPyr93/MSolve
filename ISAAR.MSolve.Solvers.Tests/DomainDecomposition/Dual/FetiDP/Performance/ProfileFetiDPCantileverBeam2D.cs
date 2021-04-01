using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Utilities;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.Performance
{
    public static class ProfileFetiDPCantileverBeam2D
    {
        public enum Precond { Dirichlet, DirichletDiagonal, Lumped }       
        public enum MatrixFormat { Dense, Skyline, SuiteSparse }

        private const double domainLengthX = 3.0, domainLengthY = 1.5;
        private const int numElementsX = 600, numElementsY = 300;
        private const int numSubdomainsX = 6, numSubdomainsY = 3;
        private const int numPcgIterations = 300;

        public static void Run()
        {
            double stiffnessRatio = 1.0;
            Precond precond = Precond.Dirichlet;
            MatrixFormat format = MatrixFormat.SuiteSparse;

            //TODO: For some reason this does not work with elements 600x300, subdomains 6x3, suitesparse, Dirichelt
            //LinearAlgebra.LibrarySettings.LinearAlgebraProviders = LinearAlgebra.LinearAlgebraProviderChoice.MKL;

            Model model = CreateModel(stiffnessRatio);
            (IVectorView ddDisplacements, ISolverLogger logger) = 
                SolveModelWithSubdomains(model, precond, format, stiffnessRatio == 1.0);
            Console.WriteLine(model.GlobalDofOrdering.NumGlobalFreeDofs);

            string path = @"C:\Users\Serafeim\Desktop\Profiling\logger.txt";
            logger.WriteToFile(path, "FETI-DP solver", true);
        }

        private static Model CreateModel(double stiffnessRatio)
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
            //double E1 = stiffnessRatio * E0;

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = domainLengthX;
            builder.DomainLengthY = domainLengthY;
            builder.NumSubdomainsX = numSubdomainsX;
            builder.NumSubdomainsY = numSubdomainsY;
            builder.NumTotalElementsX = numElementsX;
            builder.NumTotalElementsY = numElementsY;
            builder.YoungModulus = E0;
            //builder.YoungModuliOfSubdomains = new double[,] { { E1, E0, E0, E0 }, { E1, E0, E0, E0 } };
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 100.0);

            return builder.BuildModel();
        }

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(Model model,
            Precond precond, MatrixFormat format, bool homogeneous)
        {
            // Model
            model.ConnectDataStructures();

            // Corner nodes
            double meshTol = 1E-6;
            var cornerNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
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


            // Matrix format
            IFetiDPMatrixManagerFactory fetiMatrices = null;
            if (format == MatrixFormat.Dense) fetiMatrices = new FetiDPMatrixManagerFactoryDense();
            else if (format == MatrixFormat.Skyline)
            {
                IReorderingAlgorithm reordering = null;
                reordering = new OrderingAmdSuiteSparse();
                fetiMatrices = new FetiDPMatrixManagerFactorySkyline(reordering);
            }
            else if (format == MatrixFormat.SuiteSparse) fetiMatrices = new FetiDPMatrixManagerFactorySuitesparse();
            else throw new NotImplementedException();

            // Solver
            var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
            //solverBuilder.ProblemIsHomogeneous = false;

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
            //solverBuilder.PcgSettings = new PcgSettings() { ConvergenceTolerance = 1E-5 };
            solverBuilder.PcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = 1E-100,
                MaxIterationsProvider = new FixedMaxIterationsProvider(numPcgIterations)
            };

            FetiDPSolverSerial fetiSolver = solverBuilder.Build(model, cornerNodeSelection);
            //if (residualIsExact) exactResidualConvergence.FetiSolver = fetiSolver;

            // Run the analysis
            FetiDPSolverSerialTests.RunAnalysis(model, fetiSolver); //check dof separator

            // Gather the global displacements
            Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();

            return (globalDisplacements, fetiSolver.Logger);
        }

        
    }
}
