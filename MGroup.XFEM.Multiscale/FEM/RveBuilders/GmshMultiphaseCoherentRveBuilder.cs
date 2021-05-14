using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.XFEM.FEM.Input;
using MGroup.XFEM.FEM.Mesh;
using MGroup.XFEM.FEM.Mesh.GMSH;

namespace MGroup.XFEM.Multiscale.FEM.RveBuilders
{
    public class GmshMultiphaseCoherentRveBuilder : IRVEbuilder
    {
        private readonly PreprocessingMesh mesh;
        private readonly double rveVolume;
        private readonly int[] boundaryNodeIDs;
        private readonly Dictionary<int, IContinuumMaterial3D> phaseMaterials;

        private GmshMultiphaseCoherentRveBuilder(PreprocessingMesh mesh, double rveVolume, int[] boundaryNodeIDs,
            Dictionary<int, IContinuumMaterial3D> phaseMaterials)
        {
            this.mesh = mesh;
            this.rveVolume = rveVolume;
            this.boundaryNodeIDs = boundaryNodeIDs;
            this.phaseMaterials = phaseMaterials;
        }

        /// <summary>
        /// Initializes a new instance of <see cref="GmshMultiphaseCoherentRveBuilder"/>
        /// </summary>
        /// <param name="gmshMeshFilePath"></param>
        /// <param name="phaseMaterials">
        /// Keys are the phase ids. These must be identical to the physical group tags of elements in 
        /// <paramref name="gmshMeshFilePath"/>.
        /// </param>
        public static GmshMultiphaseCoherentRveBuilder CreateBuilder(double[] minCoords, double[] maxCoords,
            string gmshMeshFilePath, Dictionary<int, IContinuumMaterial3D> phaseMaterials, double meshTolerance = double.NaN)
        {
            var gmshReader = new GmshReader(gmshMeshFilePath, 3);
            PreprocessingMesh mesh = gmshReader.CreateMesh();
            if (double.IsNaN(meshTolerance))
            {
                meshTolerance = 1E-6 * Vector.CreateFromArray(maxCoords).Subtract(Vector.CreateFromArray(minCoords)).MinAbsolute();
            }
            int[] boundaryNodeIDs = FindBoundaryNodeIDs(mesh, minCoords, maxCoords, meshTolerance);
            double rveVolume = (maxCoords[0] - minCoords[0]) * (maxCoords[1] - minCoords[1]) * (maxCoords[2] - minCoords[2]);
            return new GmshMultiphaseCoherentRveBuilder(mesh, rveVolume, boundaryNodeIDs, phaseMaterials);
        }

        public IRVEbuilder Clone(int a) 
            => new GmshMultiphaseCoherentRveBuilder(mesh, rveVolume, boundaryNodeIDs, phaseMaterials);

        public Tuple<Model, Dictionary<int, Node>, double> GetModelAndBoundaryNodes()
        {
            var modelCreator = new ModelCreator();
            Model model = modelCreator.CreateModel3D(mesh, phaseMaterials);
            var boundaryNodes = new Dictionary<int, Node>(boundaryNodeIDs.Length);
            foreach (int n in boundaryNodeIDs) boundaryNodes[n] = model.NodesDictionary[n];
            return new Tuple<Model, Dictionary<int, Node>, double>(model, boundaryNodes, rveVolume);
        }

        private static int[] FindBoundaryNodeIDs(PreprocessingMesh mesh, double[] minCoords, double[] maxCoords, 
            double meshTolerance)
        {
            var boundaryNodeIDs = new List<int>();
            foreach (Vertex vertex in mesh.Vertices)
            {
                bool isBoundary = false;
                isBoundary |= Math.Abs(vertex.Coords[0] - minCoords[0]) <= meshTolerance;
                isBoundary |= Math.Abs(vertex.Coords[1] - minCoords[1]) <= meshTolerance;
                isBoundary |= Math.Abs(vertex.Coords[2] - minCoords[2]) <= meshTolerance;
                isBoundary |= Math.Abs(vertex.Coords[0] - maxCoords[0]) <= meshTolerance;
                isBoundary |= Math.Abs(vertex.Coords[1] - maxCoords[1]) <= meshTolerance;
                isBoundary |= Math.Abs(vertex.Coords[2] - maxCoords[2]) <= meshTolerance;
                if (isBoundary) boundaryNodeIDs.Add(vertex.ID); 
            }
            return boundaryNodeIDs.ToArray();
        }

    }
}
