using System;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace MGroup.Stochastic
{
    public class NeuralNetworkTrainedMaterial : IIsotropicContinuumMaterial3D
    {
        private readonly double[] strains = new double[6];
        private readonly double[] stresses = new double[6];
        private Matrix constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }
        private double[] stressesNew = new double[6];
        private double[] strainsNew = new double[6];
        public NeuralNetwork neuralNetwork = new NeuralNetwork();
        public double[] ConstParameters { get; set; }

        public NeuralNetworkTrainedMaterial()
        {
            neuralNetwork.ExtractNeuralNetworkParametersFromMatlab();
            neuralNetwork.InitializePreprocessing();
            neuralNetwork.constParameters = this.ConstParameters;
        }

        private Matrix GetConstitutiveMatrix()
        {
            return neuralNetwork.CalculateNeuralNetworkJacobian(strainsNew);
        }

        private void CalculateNextStressStrainPoint()
        {
            this.stressesNew = neuralNetwork.CalculateNeuralNetworkOutput(strainsNew).CopyToArray();
        }

        #region IFiniteElementMaterial Members

        public int ID { get; set; }

        public bool Modified => false;

        public void ResetModified() { }

        #endregion

        #region IFiniteElementMaterial3D Members

        public double[] Stresses => stressesNew;

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
                return constitutiveMatrix;
            }
        }

        public void UpdateMaterial(double[] dstrains)
        {
            //throw new NotImplementedException();
            this.strainsNew.CopyFrom(dstrains);
            for (int i = 0; i < strains.Length; i++)
            {
                this.strainsNew[i] += strains[i];
            }
            constitutiveMatrix = GetConstitutiveMatrix();
            this.CalculateNextStressStrainPoint();

        }

        public void ClearState()
        {
            //constitutiveMatrix.Clear();
            strains.Clear();
            stresses.Clear();
            strainsNew.Clear();
            stressesNew.Clear();
        }

        public void SaveState()
        {
            stresses.CopyFrom(stressesNew); 
            strains.CopyFrom(strainsNew);
        }

        public void ClearStresses()
        {
            stresses.Clear();
            stressesNew.Clear();
        }

        #endregion

        #region ICloneable Members

        object ICloneable.Clone() => Clone();

        public NeuralNetworkTrainedMaterial Clone()
        {
            return new NeuralNetworkTrainedMaterial() {ConstParameters = this.ConstParameters };
        }

        #endregion

    }

}
