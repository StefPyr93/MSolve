using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.ObjectManagers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.LinearSystems;


namespace ISAAR.MSolve.MSAnalysis.remoteMatImplementations
{
    public class RemoteMaterial : IContinuumMaterial3DDefGrad
    {
        private IMaterialManager materialManager;
        public RemoteMaterial(IMaterialManager materialManager)
        {
            this.materialManager = materialManager;
        }
        public double[] Stresses
        {
            get
            {
                if (materialManager.GetMaterialStress(this) == null)
                {
                    return new double[6];
                }
                else
                {
                    return materialManager.GetMaterialStress(this);
                }
            }
        }

        public IMatrixView ConstitutiveMatrix { get { return materialManager.GetMaterialConstitutiveMatrix(this); } }

        public int ID { get { return 1000; } }

        private bool modified;
        public bool Modified { get { return modified; } }

        public double[] Coordinates { get; set; }

        public double YoungModulus => throw new NotImplementedException();

        public double PoissonRatio => throw new NotImplementedException();

        public void ClearState() { }

        public void ClearStresses() { }
        
        public object Clone()
        {
            IContinuumMaterial3DDefGrad material = new RemoteMaterial(materialManager);
            materialManager.AddMaterial(material);
            return material;
        }

        public void ResetModified()
        {
            modified = false;
        }

        public void SaveState()
        {
            throw new NotImplementedException("material save state will be handled in an other way");
            //materialManager handles this function.
            
        }

        public void UpdateMaterial(double[] strains)
        {
            materialManager.UpdateMaterialStrainForRemote(this, strains);
        }
    }
}
