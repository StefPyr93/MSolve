using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual
{
    public static class Utilities
    {
        internal static void AnalyzeMultiSubdomainModel(IModel model, ISolverMpi solver)
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

        internal static (IMatrixView matrix, IVectorView rhs, IVectorView sol) AnalyzeSingleSubdomainModel(Model model, 
            bool suiteSparse)
        {
            int singleSubdomainID = 0;
            Utilities.RemoveSubdomains(model, singleSubdomainID);

            // Solver
            ISolver solver = null;
            if (suiteSparse)
            {
                var solverBuilder = new SuiteSparseSolver.Builder();
                solver = solverBuilder.BuildSolver(model);
                solver.PreventFromOverwrittingSystemMatrices();
            }
            else
            {
                var solverBuilder = new SkylineSolver.Builder();
                solver = solverBuilder.BuildSolver(model);
                solver.PreventFromOverwrittingSystemMatrices();
            }
            

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            ILinearSystem sys = solver.LinearSystems[singleSubdomainID];
            return  (sys.Matrix, sys.RhsVector, sys.Solution);
        }

        public static void RemoveSubdomains(Model model, int singleSubdomainID = 0)
        {
            // Replace the existing subdomains with a single one 
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.ElementsDictionary.Values) subdomain.Elements.Add(element.ID, element);
        }
    }
}
