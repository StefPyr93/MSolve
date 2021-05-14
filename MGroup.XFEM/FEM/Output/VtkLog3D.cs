using System;
using System.Linq;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using MGroup.XFEM.Output.Vtk;

namespace MGroup.XFEM.FEM.Output
{
    /// <summary>
    /// Observer that logs the displacement, strain and stress field at each analysis step to .vtk output files (1 file per 
    /// analysis step).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkLog3D : IAnalyzerLog
    {
        private readonly Model model;
        private readonly string pathNoExtension;
        private readonly bool logDisplacements, logStrains, logStresses;
        private readonly ContinuousOutputMesh vtkMesh;

        //TODO: This should be controlled by the Analyzer and passed in through IAnalyzerLog.StoreResults().
        private int iteration;

        /// <summary>
        /// Instantiates a new log that writes the displacement field of each iteration to output files, which can then be 
        /// processed in Paraview.
        /// </summary>
        /// <param name="directory">The full path of the folder whre the output files will be written, e.g. 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview.</param>
        /// <param name="filename">The name of the displacement file(s) without extensions, e.g. displacements. Keep in mind
        ///     that the actual files will be suffixed with the iteration number, e.g. 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview\\displacements_0.vtk, 
        ///     C:\\Users\\MyUser\\Desktop\\Paraview\\displacements_1.vtk, etc.</param>
        /// <param name="model">A collection of nodes, elements and other entities.</param>
        /// <param name="vtkMesh">The same mesh that is contained in <paramref name="model"/>, but expressed in VTK objects. This 
        ///     object should be shared across other VTK logs to reduce memory consumption.</param>
        /// <param name="logDisplacements">If true, the displacement field will also be written to the output file.</param>
        /// <param name="logStrains">If true, the strain field will also be written to the output file.</param>
        /// <param name="logStresses">If true, the stress field will also be written to the output file.</param>
        /// <param name="vonMisesStressCalculator">The strategy used for calculating von Mises equivalent stress. If null is 
        ///     passed, then von Mises equivalent stresses will not be written to the output file.</param>
        public VtkLog3D(string directory, string filename, Model model, ContinuousOutputMesh vtkMesh, 
            bool logDisplacements, bool logStrains, bool logStresses)
        {
            this.pathNoExtension = directory + "\\" + filename;
            this.model = model;
            this.vtkMesh = vtkMesh;
            this.logDisplacements = logDisplacements;
            this.logStrains = logStrains;
            this.logStresses = logStresses;
            iteration = 0;
        }

        public void StoreResults(DateTime startTime, DateTime endTime, IVectorView solution)
        {
            string path = pathNoExtension + $"_{iteration}.vtk";
            using (var writer = new VtkFileWriter(path))
            {
                writer.WriteMesh(vtkMesh);

                int numPoints = vtkMesh.NumOutVertices;

                // Find the displacements of each VTK point and write them to the output file
                if (logDisplacements)
                {
                    var displacementField = new DisplacementField3D(model);
                    displacementField.FindNodalDisplacements(solution);
                    var displacements = new double[numPoints][]; //TODO: this conversion between data structures should be avoided. Redundant memory and computations.
                    for (int i = 0; i < numPoints; ++i)
                    {
                        // 1-1 correspondance between nodes and VTK points, but VTK points have IDs that start from 0.
                        displacements[vtkMesh.VerticesList[i].ID] = displacementField[vtkMesh.OriginalVerticesList[i]];
                    }
                    writer.WriteVector3DField("displacements", displacements);
                }

                if (logStrains || logStresses)
                {
                    var tensorsField = new StrainStressField3D(model);
                    tensorsField.CalculateNodalTensors(solution);

                    if (logStrains)
                    {
                        var strains = new double[numPoints][]; //TODO: this conversion between data structures should be avoided. Redundant memory and computations.
                        for (int i = 0; i < numPoints; ++i)
                        {
                            strains[vtkMesh.VerticesList[i].ID] = tensorsField.GetStrainsOfNode(vtkMesh.OriginalVerticesList[i]);
                        }
                        writer.WriteTensor3DField("strain", strains);
                    }

                    if (logStresses)
                    {
                        var stresses = new double[numPoints][];
                        for (int i = 0; i < numPoints; ++i)
                        {
                            stresses[vtkMesh.VerticesList[i].ID] = tensorsField.GetStressesOfNode(vtkMesh.OriginalVerticesList[i]);
                        }
                        writer.WriteTensor3DField("stress", stresses);
                    }
                }
            }

            ++iteration;
        }
    }
}
