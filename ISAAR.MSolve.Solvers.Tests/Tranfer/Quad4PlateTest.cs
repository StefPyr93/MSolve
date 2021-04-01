using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.Transfer
{
    public static class Quad4PlateTest
    {
        //[Fact]
        private static void RunTest()
        {
            Model model = CreateModel();

            model.ConnectDataStructures();

            // Serialize each subdomain
            int numSubdomains = model.SubdomainsDictionary.Count;
            IReadOnlyList<Subdomain> originalSubdomains = model.SubdomainsDictionary.Values.ToArray();
            var serializedSubdomains = new SubdomainDto[numSubdomains];
            for (int s = 0; s < numSubdomains; ++s)
            {
                serializedSubdomains[s] = SubdomainDto.Serialize(originalSubdomains[s], model.DofSerializer);
            }

            // Deserialize each subdomain
            var deserializedSubdomains = new Subdomain[numSubdomains];
            for (int s = 0; s < numSubdomains; ++s)
            {
                deserializedSubdomains[s] = serializedSubdomains[s].Deserialize(model.DofSerializer);
            }

            // Order dofs
            foreach (Subdomain subdomain in deserializedSubdomains)
            {
                subdomain.ConnectDataStructures();
                var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
                subdomain.FreeDofOrdering = dofOrderer.OrderFreeDofs(subdomain);
            }

            // Create linear systems
            var linearSystems = new List<SingleSubdomainSystem<SkylineMatrix>>();
            foreach (Subdomain subdomain in deserializedSubdomains)
            {
                var ls = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
                ls.Reset();
                subdomain.Forces = Vector.CreateZero(ls.Size);
                linearSystems.Add(ls);
            }


            // Create the stiffness matrices
            foreach (ILinearSystem ls in linearSystems)
            {
                ISubdomain subdomain = ls.Subdomain;
                var provider = new ElementStructuralStiffnessProvider();
                var assembler = new SkylineAssembler();
                SkylineMatrix stiffness = assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, 
                    subdomain.EnumerateElements(), provider);
                ls.Matrix = stiffness;
            }


            // Calculate the trace of each stiffness matrix
            double[] tracesComputed = linearSystems.Select(ls => ls.Matrix.Trace()).ToArray();

            // Check correctness
            // double[] tracesExpected = { 1.30846153846154E-05, 1.64855769230769E-05, 1.99384615384615E-05, 2.32615384615385E-05 };
            double[] tracesExpected = CalcReferenceStiffnessTraces();
            double tol = 1E-13;
            Assert.True(Vector.CreateFromArray(tracesExpected).Equals(Vector.CreateFromArray(tracesComputed), tol));

        }

        public static double[] CalcReferenceStiffnessTraces()
        {
            Model model = CreateModel();

            model.ConnectDataStructures();

            // Order dofs
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            dofOrderer.OrderFreeDofs(model);

            // Create linear systems
            var linearSystems = new List<SingleSubdomainSystem<SkylineMatrix>>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                var ls = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
                linearSystems.Add(ls);
                ls.Reset();
                subdomain.Forces = Vector.CreateZero(ls.Size);
            }

            // Create the stiffness matrices
            var provider = new ElementStructuralStiffnessProvider();
            foreach (ILinearSystem ls in linearSystems)
            {
                ISubdomain subdomain = ls.Subdomain;
                var assembler = new SkylineAssembler();
                SkylineMatrix stiffness = assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, 
                    subdomain.EnumerateElements(), provider);
                ls.Matrix = stiffness;
            }

            // Calculate the trace of each stiffness matrix
            double[] traces = linearSystems.Select(ls => ls.Matrix.Trace()).ToArray();
            return traces;
        }

        public static Model CreateModel()
        {
            double E = 2.1E-7;
            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 2.0;
            builder.DomainLengthY = 2.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 8;
            builder.NumTotalElementsY = 8;
            builder.YoungModuliOfSubdomains = new double[,] { { E, 1.25 * E }, { 1.5 * E, 1.75 * E } };
            //builder.YoungModulus = E;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, StructuralDof.TranslationY, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerRightCorner, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationX, 100.0);

            return builder.BuildModel();
        }

        public static double Trace(this IMatrixView matrix)
        {
            double trace = 0.0;
            for (int i = 0; i < matrix.NumRows; ++i) trace += matrix[i, i];
            return trace;
        }
    }
}
