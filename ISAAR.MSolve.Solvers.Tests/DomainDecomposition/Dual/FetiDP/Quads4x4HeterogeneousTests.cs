using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4HeterogeneousTests
    {
        [Fact]
        public static void TestScalingMatrices() //TODO: Most of the code is the same as the homogeneous case. Remove duplication.
        {
            #region Replace the next with hardcoded matrices and mocking objects
            // Run the analysis so that all objects are created
            // Setup the model
            double stiffnessRatio = 1E-2; // Do not change this! The expected solution is taken for this value
            Model model = Example4x4QuadsHeterogeneous.CreateModel(stiffnessRatio);

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolverOLD.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHeterogeneous.DefineCornerNodeSelectionSerial(model);
            var fetiSolverBuilder = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            var preconditionerFactory = new DirichletPreconditionerOLD.Factory();
            fetiSolverBuilder.PreconditionerFactory = preconditionerFactory;
            FetiDPSolverOLD fetiSolver = fetiSolverBuilder.BuildSolver(model);

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();
            #endregion

            // Access private fields of FetiDPSolver and DirichletPreconditioner.Factory using reflection
            FieldInfo fi = typeof(FetiDPSolverOLD).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumeratorOLD)fi.GetValue(fetiSolver);
            fi = typeof(FetiDPSolverOLD).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparatorOLD)fi.GetValue(fetiSolver);
            fi = typeof(FetiDPSolverOLD).GetField("stiffnessDistribution", BindingFlags.NonPublic | BindingFlags.Instance);
            var stiffnessDistribution = (IStiffnessDistributionOLD)fi.GetValue(fetiSolver);
            MethodInfo method = preconditionerFactory.GetType().GetMethod("CalcBoundaryPreconditioningBooleanMatrices",
                BindingFlags.NonPublic | BindingFlags.Instance);
            var Bpbr = (Dictionary<int, IMappingMatrix>)method.Invoke(preconditionerFactory, 
                new object[] { model, stiffnessDistribution, dofSeparator, lagrangeEnumerator });

            // Compare the mapping matrices against the expected ones
            double tol = 1E-13;
            for (int s = 0; s < 4; ++s)
            {
                Matrix explicitBpr = Bpbr[s].MultiplyRight(Matrix.CreateIdentity(Bpbr[s].NumColumns));
                Assert.True(Example4x4QuadsHeterogeneous.GetMatrixBpbr(s).Equals(explicitBpr, tol));
            }
        }

        [Fact]
        public static void TestSolver()
        {
            // Setup the model
            double stiffnessRatio = 1E-2; // Do not change this! The expected solution is taken for this value
            Model model = Example4x4QuadsHeterogeneous.CreateModel(stiffnessRatio);

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolverOLD.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHeterogeneous.DefineCornerNodeSelectionSerial(model);
            var fetiSolverBuilder = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditionerOLD.Factory();
            FetiDPSolverOLD fetiSolver = fetiSolverBuilder.BuildSolver(model);

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);

            // Check against expected solution
            double tol = 1E-7;
            Assert.True(Example4x4QuadsHeterogeneous.SolutionGlobalDisplacements.Equals(globalU, tol));
        }
    }
}
