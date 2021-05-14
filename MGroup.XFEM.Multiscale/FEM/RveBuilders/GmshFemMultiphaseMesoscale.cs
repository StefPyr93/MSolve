using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;

namespace MGroup.XFEM.Multiscale.FEM.RveBuilders
{
    public class GmshFemMultiphaseMesoscale : IMesoscale
    {
        private const int dimension = 3;

        #region input 
        public double[] CoordsMin { get; set; }
        public double[] CoordsMax { get; set; }

        public IContinuumMaterial3D MatrixMaterial { get; set; }

        public IContinuumMaterial3D InclusionMaterial { get; set; }

        public ISolverBuilder SolverBuilder { get; set; } = new SuiteSparseSolver.Builder();

        public string GmshMeshFilePath { get; set; }

        public double[] TotalStrain { get; set; }

        public int NumLoadingIncrements { get; set; } = 10;
        #endregion

        #region output
        public IList<double[]> Strains { get; set; } = new List<double[]>();

        public IList<double[]> Stresses { get; set; } = new List<double[]>();

        public IList<IMatrixView> ConstitutiveMatrices { get; set; } = new List<IMatrixView>();
        #endregion

        public void RunAnalysis()
        {
            var phaseMaterials = new Dictionary<int, IContinuumMaterial3D>();
            phaseMaterials[1] = MatrixMaterial;
            phaseMaterials[2] = InclusionMaterial;
            var rveBuilder = GmshMultiphaseCoherentRveBuilder.CreateBuilder(
                CoordsMin, CoordsMax, GmshMeshFilePath, phaseMaterials);
            var microstructure = new Microstructure3D(rveBuilder, model => SolverBuilder.BuildSolver(model), false, 1);

            for (int i = 0; i < NumLoadingIncrements; ++i)
            {
                double[] macroStrain = TotalStrain.Scale((i + 1) / (double)NumLoadingIncrements);
                microstructure.UpdateMaterial(macroStrain);

                Strains.Add(macroStrain.Copy());
                Stresses.Add(microstructure.Stresses.Copy());
                ConstitutiveMatrices.Add(microstructure.ConstitutiveMatrix.Copy());

                microstructure.SaveState();
            }
        }

    }
}
