using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Analyzers.ObjectManagers
{
    public interface IMaterialManager
    {
        void Initialize();

        void AddMaterial(IContinuumMaterial3DDefGrad addedMaterial);
        void UpdateMaterialStrainForRemote(IContinuumMaterial3DDefGrad ghostMaterial, double[] ghostMaterialStrain);
        double[] GetMaterialStress(IContinuumMaterial3DDefGrad ghostMaterial);
        IMatrixView GetMaterialConstitutiveMatrix(IContinuumMaterial3DDefGrad ghostMaterial);
        //private void BuildMaterials();
        void UpdateMaterials();
        void SaveState();
        
    }
}
