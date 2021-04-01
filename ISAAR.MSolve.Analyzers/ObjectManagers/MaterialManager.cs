using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using System.Linq;

namespace ISAAR.MSolve.Analyzers.ObjectManagers
{
    public class MaterialManager: IMaterialManager
    {
        private IContinuumMaterial3DDefGrad chosenMaterial;

        private List<IContinuumMaterial3DDefGrad> pseudoMaterials = new List<IContinuumMaterial3DDefGrad>();

        private Dictionary<IContinuumMaterial3DDefGrad, double[]> pseudoMaterialStrains;
        private Dictionary<IContinuumMaterial3DDefGrad, double[]> pseudoMaterialStresses;
        private Dictionary<IContinuumMaterial3DDefGrad, IMatrixView> pseudoMaterialConsMatrices;

        private Dictionary<IContinuumMaterial3DDefGrad, int> pseudoMaterialsMappingToDatabase;

        private IContinuumMaterial3DDefGrad[] materialDatabase;

        public MaterialManager(IContinuumMaterial3DDefGrad coosenMaterial)
        {
            this.chosenMaterial = coosenMaterial;
        }

        public void AddMaterial(IContinuumMaterial3DDefGrad addedMaterial)
        {
            pseudoMaterials.Add(addedMaterial);
        }

        public void UpdateMaterialStrainForRemote(IContinuumMaterial3DDefGrad pseudoMaterial, double[] pseudoMaterialStrain)
        {
            pseudoMaterialStrains[pseudoMaterial] = pseudoMaterialStrain;
        }

        public double[] GetMaterialStress(IContinuumMaterial3DDefGrad pseudoMaterial)
        {
            return pseudoMaterialStresses[pseudoMaterial];
        }

        public IMatrixView GetMaterialConstitutiveMatrix(IContinuumMaterial3DDefGrad pseudoMaterial)
        {
            return pseudoMaterialConsMatrices[pseudoMaterial];
        }
        public void Initialize()
        {
            BuildMaterials();

            pseudoMaterialStrains = pseudoMaterials.Select(x => new KeyValuePair<IContinuumMaterial3DDefGrad, double[]>(x, null)).ToDictionary(x =>x.Key,x=>x.Value);
            pseudoMaterialStresses = pseudoMaterials.Select(x => new KeyValuePair<IContinuumMaterial3DDefGrad, double[]>(x, null)).ToDictionary(x => x.Key, x => x.Value);
            pseudoMaterialConsMatrices = pseudoMaterials.Select(x=>  new KeyValuePair<IContinuumMaterial3DDefGrad,IMatrixView>( x,null)).ToDictionary(x=>x.Key,x=>x.Value);


            foreach(var ghoMat in pseudoMaterials)
            {
                pseudoMaterialConsMatrices[ghoMat] = materialDatabase[pseudoMaterialsMappingToDatabase[ghoMat]].ConstitutiveMatrix; //.CopyToFullMatrix();
            }

            //...
        }

        private void BuildMaterials()
        {
            materialDatabase = new IContinuumMaterial3DDefGrad[pseudoMaterials.Count];
            pseudoMaterialsMappingToDatabase = new Dictionary<IContinuumMaterial3DDefGrad, int>(pseudoMaterials.Count);
            int counter = 0;
            foreach(var material in pseudoMaterials)
            {
                materialDatabase[counter] = (IContinuumMaterial3DDefGrad)chosenMaterial.Clone();
                pseudoMaterialsMappingToDatabase[material] = counter; counter++;
            }
        }

        public void UpdateMaterials()
        {
            foreach(KeyValuePair<IContinuumMaterial3DDefGrad,double[]> matAndStrain  in pseudoMaterialStrains)
            {
                materialDatabase[pseudoMaterialsMappingToDatabase[matAndStrain.Key]].UpdateMaterial(matAndStrain.Value);
            }
            foreach (var ghoMat in pseudoMaterials)
            {
                pseudoMaterialConsMatrices[ghoMat] = materialDatabase[pseudoMaterialsMappingToDatabase[ghoMat]].ConstitutiveMatrix; //.CopyToFullMatrix();
            }
            foreach (var ghoMat in pseudoMaterials)
            {
                pseudoMaterialStresses[ghoMat] = materialDatabase[pseudoMaterialsMappingToDatabase[ghoMat]].Stresses; //.CopyToFullMatrix();
            }
        }

        public void SaveState()
        {
            foreach (KeyValuePair<IContinuumMaterial3DDefGrad, double[]> matAndStrain in pseudoMaterialStrains)
            {
                materialDatabase[pseudoMaterialsMappingToDatabase[matAndStrain.Key]].SaveState();
            }
        }


    }
}
