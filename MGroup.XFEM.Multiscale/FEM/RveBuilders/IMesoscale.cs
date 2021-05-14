using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Solvers;

namespace MGroup.XFEM.Multiscale.FEM.RveBuilders
{
    public interface IMesoscale
    {
        #region input 
        //int Seed { get; set; }

        //double VolumeFraction { get; set; }

        IContinuumMaterial3D MatrixMaterial { get; set; }

        IContinuumMaterial3D InclusionMaterial { get; set; }

        ISolverBuilder SolverBuilder { get; set; }

        double[] TotalStrain { get; set; }

        int NumLoadingIncrements { get; set; }
        #endregion

        #region output
        IList<double[]> Strains { get; set; }

        IList<double[]> Stresses { get; set; }

        IList<IMatrixView> ConstitutiveMatrices { get; set; }
        #endregion

        void RunAnalysis();
    }
}
