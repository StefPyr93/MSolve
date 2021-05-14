using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MGroup.XFEM.FEM.Elements;
using MGroup.XFEM.FEM.Input;
using MGroup.XFEM.FEM.Mesh;
using MGroup.XFEM.FEM.Mesh.GMSH;

namespace MGroup.XFEM.Tests.FEM
{
    public static class FemUtilities
    {
        public static void ApplyBCsCantileverTension(Model model, int dim)
        {
            // Boundary conditions
            double meshTol = 1E-7;

            // Left side: Ux=Uy=Uz=0
            double minX = model.NodesDictionary.Values.Select(n => n.X).Min();
            foreach (var node in model.NodesDictionary.Values.Where(n => Math.Abs(n.X - minX) <= meshTol))
            {
                node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0 });
                if (dim >= 2) node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0 });
                if (dim == 3) node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ, Amount = 0 });
            }

            // Right side: Fx = 100
            double maxX = model.NodesDictionary.Values.Select(n => n.X).Max();
            Node[] rightSideNodes = model.NodesDictionary.Values.Where(n => Math.Abs(n.X - maxX) <= meshTol).ToArray();
            double load = 1.0 / rightSideNodes.Length;
            foreach (var node in rightSideNodes)
            {
                model.Loads.Add(new Load() { Node = node, DOF = StructuralDof.TranslationX, Amount = load });
            }
        }

        public static Model Create3DModelFromGmsh(string gmshMeshFile, 
            Dictionary<int, IContinuumMaterial3D> materialsOfPhysicalGroups)
        {
            var gmshReader = new GmshReader(gmshMeshFile, 3);
            PreprocessingMesh mesh = gmshReader.CreateMesh();
            var modelCreator = new ModelCreator();
            return modelCreator.CreateModel3D(mesh, materialsOfPhysicalGroups);
        }

        public static IVectorView RunStaticLinearAnalysis(Model model, ILogFactory logFactory = null, 
            ISolverBuilder solverBuilder = null)
        {
            Console.WriteLine("Starting analysis");
            if (solverBuilder == null) solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);
            var problem = new ProblemStructural(model, solver);
            var linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);

            // Output
            if (logFactory != null) linearAnalyzer.LogFactories[0] = logFactory;
            

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            Console.WriteLine("Analysis finished");
            return solver.LinearSystems[0].Solution;
        }
    }
}
