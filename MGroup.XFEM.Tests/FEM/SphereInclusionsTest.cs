using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using MGroup.XFEM.FEM.Output;
using MGroup.XFEM.Multiscale.FEM.RveBuilders;
using Xunit;

namespace MGroup.XFEM.Tests.FEM
{
    /// <summary>
    /// For the .geo script that creates the mesh, see Resources\gmsh_sphere_inclusions\spheres_in_rve.geo
    /// </summary>
    public static class SphereInclusionsTest
    {
        private const string workingDirectory = @"C:\Users\stefp\OneDrive\Desktop";
        private const string meshFile = workingDirectory + "\\spheres_in_rve.msh";

        private const double matrixE = 1E0, inclusionE = 1E3;

        [Fact]
        public static void RunLinearAnalysis()
        {
            var materialsOfPhysicalGroups = new Dictionary<int, IContinuumMaterial3D>();
            materialsOfPhysicalGroups[1] = new ElasticMaterial3D() { YoungModulus = matrixE, PoissonRatio = 0.3 };
            materialsOfPhysicalGroups[2] = new ElasticMaterial3D() { YoungModulus = inclusionE, PoissonRatio = 0.3 };
            Model model = FemUtilities.Create3DModelFromGmsh(meshFile, materialsOfPhysicalGroups);
            FemUtilities.ApplyBCsCantileverTension(model, 3);

            string outputDirectory = workingDirectory;
            ILogFactory logs = new VtkLogFactory(3, model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            FemUtilities.RunStaticLinearAnalysis(model, logFactory: logs);
        }

        [Fact]
        public static void RunMultiscaleAnalysis()
        {
            var mesoscale = new GmshFemMultiphaseMesoscale();
            mesoscale.CoordsMin = new double[] { -1, -1, -1 };
            mesoscale.CoordsMax = new double[] { +1, +1, +1 };
            mesoscale.GmshMeshFilePath = meshFile;
            mesoscale.MatrixMaterial = new ElasticMaterial3D() { YoungModulus = matrixE, PoissonRatio = 0.3 };
            mesoscale.InclusionMaterial = new ElasticMaterial3D() { YoungModulus = inclusionE, PoissonRatio = 0.3 };
            mesoscale.TotalStrain = new double[] { 0.04, 0.02, 0.02, 0.04, 0.02, 0.02 };
            mesoscale.NumLoadingIncrements = 10;

            mesoscale.RunAnalysis();
           
        }

    }
}
