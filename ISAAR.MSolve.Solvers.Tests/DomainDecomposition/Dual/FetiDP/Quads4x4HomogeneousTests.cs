using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MPI;
using Xunit;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.UnitTests;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using ISAAR.MSolve.Solvers.Logging;

//TODO: Also test stiffness distribution and preconditioners in other classes.
//TODO: Create the dofSeparator and lagrangeEnumerator manually, without using FetiDPSolver.
//TODO: TestInterfaceProblemSolution should mock matrices and vectors from TestInterfaceProblemCreation.
//TODO: There is a lot of code duplication between the methods.
//TODO: The results of matrix slicing and factorization should be hardcoded and tested against instead.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4HomogeneousTests
    {
        [Fact]
        public static void TestCoarseProblem()
        {
            // Setup the model and solver
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            var solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);
            solver.Initialize();

            // Use the hardcoded intermediate matrices & vectors
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers = MockStiffnesses();
            var fbc = new Dictionary<int, Vector>();
            var fr = new Dictionary<int, Vector>();
            for (int s = 0; s < 4; ++s)
            {
                fbc[s] = Example4x4QuadsHomogeneous.GetVectorFbc(s);
                fr[s] = Example4x4QuadsHomogeneous.GetVectorFr(s);
            }

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolverOLD).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparatorOLD)fi.GetValue(solver);

            // Calculate the coarse problem matrix and rhs
            var coarseSolver = new FetiDPCoarseProblemSolverDenseOLD(model);
            Vector globalFcStar = coarseSolver.CreateCoarseProblemRhs(dofSeparator, matrixManagers, fr, fbc);
            MethodInfo method = coarseSolver.GetType().GetMethod("CreateGlobalKccStar",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Matrix globalKccStar = (Matrix)method.Invoke(coarseSolver, new object[] { dofSeparator, matrixManagers });

            // Check against expected matrices
            var expectedKccStar = Example4x4QuadsHomogeneous.MatrixGlobalKccStar;
            var expectedFcStar = Example4x4QuadsHomogeneous.VectorGlobalFcStar;
            double tol = 1E-13;
            Assert.True(expectedKccStar.Equals(globalKccStar, tol));
            Assert.True(expectedFcStar.Equals(globalFcStar, tol));
        }

        [Fact]
        public static void TestDisconnectedDisplacements()
        {
            //TODO: Perhaps use the Br, Bc from the class that tests them instead of the solver.

            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            FetiDPSolverOLD solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);
            solver.Initialize();

            // Mock the stiffness matrices and force vectors
            var fr = new Dictionary<int, Vector>();
            for (int s = 0; s < 4; ++s)
            {
                fr[s] = Example4x4QuadsHomogeneous.GetVectorFr(s);
            }
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers = MockStiffnesses();

            // Use reflection to set the necessary matrices
            FieldInfo fi;
            fi = typeof(FetiDPSolverOLD).GetField("matrixManagers", BindingFlags.NonPublic | BindingFlags.Instance);
            fi.SetValue(solver, matrixManagers);

            Vector dr = solver.CalcDisconnectedDisplacements(fr);
            var expectedDr = Example4x4QuadsHomogeneous.VectorDr;

            double tol = 1E-13;
            Assert.True(expectedDr.Equals(dr, tol));
        }

        [Fact]
        public static void TestFlexibilityMatrices()
        {
            // Setup the model and solver
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            var solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);
            solver.Initialize();

            // Mock the stiffness matrices
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers = MockStiffnesses();

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolverOLD).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumeratorOLD)fi.GetValue(solver);
            fi = typeof(FetiDPSolverOLD).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparatorOLD)fi.GetValue(solver);

            // Create the flexibility matrices by multiplying with identity matrices
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            var flexibility = new FetiDPFlexibilityMatrixOLD(dofSeparator, lagrangeEnumerator, matrixManagers);
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numLagranges, flexibility.MultiplyFIrr);
            Matrix FIrc = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numLagranges, 
                (x, y) => y.CopyFrom(flexibility.MultiplyFIrc(x)));

            // Check against expected matrices
            double tol = 1E-11;
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrr.Equals(FIrr, tol));
            Assert.True(Example4x4QuadsHomogeneous.MatrixFIrc.Equals(FIrc, tol));
        }

        [Fact]
        public static void TestInterfaceProblemCreation()
        {
            // Setup the model and solver
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            var solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);
            solver.Initialize();

            // Mock the stiffness matrices and force vectors
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers = MockStiffnesses();
            Vector dr = Example4x4QuadsHomogeneous.VectorDr;

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolverOLD).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumeratorOLD)fi.GetValue(solver);
            fi = typeof(FetiDPSolverOLD).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparatorOLD)fi.GetValue(solver);

            // Hardcoded coarse problem matrix and rhs
            var coarseSolver = new FetiDPCoarseProblemSolverDenseOLD(model);
            Vector globalFcStar = Example4x4QuadsHomogeneous.VectorGlobalFcStar;
            Matrix inverseKccStar = Example4x4QuadsHomogeneous.MatrixGlobalKccStar.Invert(); // It must be set as a private field using reflection.
            fi = typeof(FetiDPCoarseProblemSolverDenseOLD).GetField("inverseGlobalKccStar",
                BindingFlags.NonPublic | BindingFlags.Instance);
            fi.SetValue(coarseSolver, inverseKccStar);

            // Create the rhs vector of the interface problem 
            FetiDPInterfaceProblemSolverOLD interfaceSolver = new FetiDPInterfaceProblemSolverOLD.Builder().Build();
            var flexibility = new FetiDPFlexibilityMatrixOLD(dofSeparator, lagrangeEnumerator, matrixManagers);
            Vector fcStar = Example4x4QuadsHomogeneous.VectorGlobalFcStar;
            MethodInfo method = interfaceSolver.GetType().GetMethod("CreateInterfaceProblemRhs",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector interfaceRhs = (Vector)method.Invoke(interfaceSolver, new object[] { flexibility, coarseSolver, fcStar, dr });

            // Create the matrix of the interface problem by multiplying with identity matrix
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var interfaceMatrixImplicit = 
                new FetiDPInterfaceProblemSolverOLD.InterfaceProblemMatrix(flexibility, coarseSolver);
            Matrix interfaceMatrix = ImplicitMatrixUtilities.MultiplyWithIdentity(numLagranges, numLagranges, 
                interfaceMatrixImplicit.Multiply); // Action<T> is contravariant!!!

            // Check against expected linear system
            double tol = 1E-13;
            Assert.True(Example4x4QuadsHomogeneous.InterfaceProblemRhs.Equals(interfaceRhs, tol));
            Assert.True(Example4x4QuadsHomogeneous.InterfaceProblemMatrix.Equals(interfaceMatrix, tol));
        }

        [Fact]
        public static void TestInterfaceProblemSolution()
        {
            // Setup the model and solver
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            var solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);
            solver.Initialize();

            // Mock the stiffness matrices and force vectors
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers = MockStiffnesses();
            Dictionary<int, IFetiSubdomainMatrixManagerOLD> matrixManagersPreconditioning = MockStiffnessesForPreconditioning();
            //Dictionary<int, Matrix> Krr = MatricesKrr;
            Vector dr = Example4x4QuadsHomogeneous.VectorDr;

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolverOLD).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumeratorOLD)fi.GetValue(solver);
            fi = typeof(FetiDPSolverOLD).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparatorOLD)fi.GetValue(solver);

            // Hardcoded coarse problem matrix and rhs
            var coarseSolver = new FetiDPCoarseProblemSolverDenseOLD(model);
            Vector globalFcStar = Example4x4QuadsHomogeneous.VectorGlobalFcStar;
            Matrix inverseKccStar = Example4x4QuadsHomogeneous.MatrixGlobalKccStar.Invert(); // It must be set as a private field using reflection.
            fi = typeof(FetiDPCoarseProblemSolverDenseOLD).GetField("inverseGlobalKccStar",
                BindingFlags.NonPublic | BindingFlags.Instance);
            fi.SetValue(coarseSolver, inverseKccStar);

            // Dirichlet preconditioner
            var precondFactory = new DirichletPreconditionerOLD.Factory();
            //var repackagedKrr = new Dictionary<int, IMatrixView>();
            //foreach (var idMatrixPair in Krr) repackagedKrr[idMatrixPair.Key] = idMatrixPair.Value;
            var stiffnessDistribution = new FetiDPHomogeneousStiffnessDistributionOLD(model, dofSeparator);
            stiffnessDistribution.Update(null);
            IFetiPreconditioner preconditioner = precondFactory.CreatePreconditioner(model,
                stiffnessDistribution, dofSeparator, lagrangeEnumerator, matrixManagersPreconditioning);

            // Solve the interface problem
            FetiDPInterfaceProblemSolverOLD interfaceSolver = new FetiDPInterfaceProblemSolverOLD.Builder().Build();
            var flexibility = new FetiDPFlexibilityMatrixOLD(dofSeparator, lagrangeEnumerator, matrixManagers);
            var logger = new SolverLoggerOLD("mock FETI-DP");
            (Vector lagranges, Vector uc) = interfaceSolver.SolveInterfaceProblem(
                flexibility, preconditioner, coarseSolver, globalFcStar, dr, Example4x4QuadsHomogeneous.GlobalForcesNorm, logger);

            // Check against expected solution
            double tol = 1E-7;
            Assert.True(Example4x4QuadsHomogeneous.SolutionLagrangeMultipliers.Equals(lagranges, tol));
            Assert.True(Example4x4QuadsHomogeneous.SolutionCornerDisplacements.Equals(uc, tol));
        }

        [Fact]
        public static void TestScalingMatrices()
        {
            #region Replace the next with hardcoded matrices and mocking objects
            // Run the analysis so that all objects are created
            // Setup the model
            Model model = Example4x4QuadsHomogeneous.CreateModel();

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolverOLD.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
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
                Assert.True(Example4x4QuadsHomogeneous.GetMatrixBpbr(s).Equals(explicitBpr, tol));
            }
        }

        [Fact]
        public static void TestSolver()
        {
            // Setup the model and solver
            Model model = Example4x4QuadsHomogeneous.CreateModel();
            var fetiMatrices = new FetiDPSubdomainMatrixManagerDenseOLD.Factory();
            var cornerNodeSelection = Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);
            var solver = new FetiDPSolverOLD.Builder(cornerNodeSelection, fetiMatrices).BuildSolver(model);
            var problem = new ProblemStructural(model, solver);
            var linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);

            // Run the analysis
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in solver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = solver.GatherGlobalDisplacements(sudomainDisplacements);

            // Check against expected solution
            double tol = 1E-7;
            Assert.True(Example4x4QuadsHomogeneous.SolutionGlobalDisplacements.Equals(globalU, tol));
        }

        /// <summary>
        /// Uses reflection to set the necessary matrices
        /// </summary>
        private static Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> MockStiffnesses()
        {
            var matrixManagers = new Dictionary<int, IFetiDPSubdomainMatrixManagerOLD>();
            for (int s = 0; s < 4; ++s)
            {
                var matrixManager = new FetiDPSubdomainMatrixManagerDenseOLD(null);
                matrixManagers[s] = matrixManager;

                FieldInfo fi;
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("Kcc", BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKcc(s));
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("Krc", BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKrc(s));
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("Krr", BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKrr(s));
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("inverseKrr", 
                    BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKrr(s).FactorCholesky(false));

                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("Kbb", BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKbb(s));
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("Kbi", BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, Example4x4QuadsHomogeneous.GetMatrixKbi(s));
                Matrix inverseKii = Example4x4QuadsHomogeneous.GetMatrixKii(s);
                inverseKii.InvertInPlace();
                fi = typeof(FetiDPSubdomainMatrixManagerDenseOLD).GetField("inverseKii",
                    BindingFlags.NonPublic | BindingFlags.Instance);
                fi.SetValue(matrixManager, inverseKii);
            }
            return matrixManagers;
        }

        /// <summary>
        /// Repackages the result of <see cref="MockStiffnesses"/>
        /// </summary>
        private static Dictionary<int, IFetiSubdomainMatrixManagerOLD> MockStiffnessesForPreconditioning()
        {
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> managersConcrete = MockStiffnesses();
            var managersGeneral = new Dictionary<int, IFetiSubdomainMatrixManagerOLD>();
            foreach (int s in managersConcrete.Keys) managersGeneral[s] = managersConcrete[s];
            return managersGeneral;
        }
    }
}
