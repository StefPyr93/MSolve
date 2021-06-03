using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using Accord.Math.Decompositions;
using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Stochastic
{
    public class ExponentialDamageMaterial : IIsotropicContinuumMaterial3D
    {
        private bool modified;

        public double youngModulus { get; set; }
        public double poissonRatio { get; set; }
        public double Strain_0 { get; set; }
        public double A { get; set; }
        public double B { get; set; }
        public double Veta { get; set; }

        object ICloneable.Clone() => Clone();
        //IContinuumMaterial3D IContinuumMaterial3D.Clone() => Clone();
        public ExponentialDamageMaterial Clone()
        {
            return new ExponentialDamageMaterial()
            {
                youngModulus = this.youngModulus,
                poissonRatio = this.poissonRatio,
                Strain_0 = this.Strain_0,
                A = this.A,
                B = this.B,
                Veta = this.Veta,
            };
        }

        private Matrix constitutiveMatrix;
        private Matrix inverseConstitutiveMatrix;
        private double[] strain;
        private double[] strain_prev;
        private double[] stress_eff;
        public double dmg;
        private double dmg_prev;
        private Matrix D_tan;
        private double[,] D_tan_prev;
        private double strain_eq;
        private double strain_eq_prev;
        private double[,] D_tan_f;

        private double[] stress = new double[6];

        private bool matrices_not_initialized = true;
        public void InitializeMatrices()
        {
            D_tan = Matrix.CreateZero(6, 6);
            D_tan_prev = new double[6, 6];
            D_tan_f = new double[6, 6];
            strain_eq_prev = Strain_0;
            strain = new double[6];
            strain_prev = new double[6];
            //stress = new double[6];
            stress_eff = new double[6];
            matrices_not_initialized = false;
            constitutiveMatrix = GetConstitutiveMatrix();
            inverseConstitutiveMatrix = GetInverseConstitutiveMatrix();
        }

        public void UpdateMaterial(double[] dstrain)
        {
            if (matrices_not_initialized)
            { this.InitializeMatrices(); }
            // Store previous Tangent moduli for ifmaterialsModified check
            for (int k = 0; k < 6; k++)
            {
                for (int j = 0; j < 6; j++)
                {
                    D_tan_prev[k, j] = D_tan[k, j];
                }
            }
            for (int i = 0; i < 6; i++)
                strain[i] = strain_prev[i] + dstrain[i];
            for (int i = 0; i < 6; i++)
            {
                stress_eff[i] = 0;
                for (int j = 0; j < 6; j++)
                {
                    stress_eff[i] += this.constitutiveMatrix[i, j] * strain[j];
                }
            }
            var strain_eq2 = 0.0;
            for (int i = 0; i < 6; i++)
            {
                //strain_eq2 += Math.Pow(strain[i], 2);
                strain_eq2 += 1 / youngModulus * strain[i] * stress_eff[i];
            }
            strain_eq = Math.Sqrt(strain_eq2);
            strain_eq2 += 1e-10;
            if (strain_eq <= strain_eq_prev)//&& dmg == dmg_prev
            {
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        D_tan[i, j] = (1 - dmg_prev) * this.constitutiveMatrix[i, j];
                    }
                }
                for (int i = 0; i < 6; i++)
                {
                    stress[i] = 0;
                    for (int j = 0; j < 6; j++)
                    {
                        stress[i] += D_tan[i, j] * strain[j];
                    }
                }
            }
            else//if (strain_eq > strain_eq_prev)
            {
                var dmg_temp = 1 - (1 - A) * Strain_0 / strain_eq - A * Math.Exp(-B * (strain_eq - Strain_0));
                if (dmg_temp > dmg_prev)
                {
                    dmg = dmg_temp;
                }

                if (dmg < 0)
                    dmg = 0;
                else if (dmg > 0.99)
                    dmg = 0.99;
                var Dtan_coeff = (1 - A) * Strain_0 / Math.Pow(strain_eq, 2) + A * B * Math.Exp(-B * (strain_eq - Strain_0));
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        //D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j] - (Dtan_coeff * stress_eff[i]) * (strain[j] / strain_eq);
                        D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j] - Dtan_coeff / (youngModulus * strain_eq) * stress_eff[i] * stress_eff[j] ;
                        D_tan_f[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j];
                        //D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j];
                    }
                }
                for (int i = 0; i < 6; i++)
                {
                    stress[i] = 0;
                    for (int j = 0; j < 6; j++)
                    {
                        stress[i] += D_tan_f[i, j] * strain[j];
                        //stress[i] += D_tan[i, j] * strain[j];
                    }
                }
            }
            this.modified = CheckIfConstitutiveMatrixChanged();
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (Math.Abs(D_tan_prev[i, j] - D_tan[i, j]) > 1e-10)
                        return true;

            return false;
        }

        private Matrix GetConstitutiveMatrix()
        {
            double fE1 = YoungModulus / (double)(1 + PoissonRatio);
            double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
            double fE3 = fE1 + fE2;
            double fE4 = fE1 * 0.5;
            var afE = Matrix.CreateZero(6, 6);
            afE[0, 0] = fE3;
            afE[0, 1] = fE2;
            afE[0, 2] = fE2;
            afE[1, 0] = fE2;
            afE[1, 1] = fE3;
            afE[1, 2] = fE2;
            afE[2, 0] = fE2;
            afE[2, 1] = fE2;
            afE[2, 2] = fE3;
            afE[3, 3] = fE4;
            afE[4, 4] = fE4;
            afE[5, 5] = fE4;

            return afE;
        }

        private Matrix GetInverseConstitutiveMatrix()
        {
            double fE1 = 1 / YoungModulus;
            double fE2 = fE1 * -PoissonRatio;
            double fE3 = fE1 * 2 * (1 + PoissonRatio);
            var afE = Matrix.CreateZero(6, 6);
            afE[0, 0] = fE1;
            afE[0, 1] = fE2;
            afE[0, 2] = fE2;
            afE[1, 0] = fE2;
            afE[1, 1] = fE1;
            afE[1, 2] = fE2;
            afE[2, 0] = fE2;
            afE[2, 1] = fE2;
            afE[2, 2] = fE1;
            afE[3, 3] = fE3;
            afE[4, 4] = fE3;
            afE[5, 5] = fE3;

            return afE;
        }

        public double[] Stresses
        {
            get { return stress; }
        }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (D_tan == null) UpdateMaterial(new double[6]);
                return D_tan;
            }
        }

        public void SaveState()
        {
            if (strain_eq > strain_eq_prev)
                strain_eq_prev = strain_eq;
            dmg_prev = dmg;
            strain_prev = new double[6] { strain[0], strain[1], strain[2], strain[3], strain[4], strain[5] };
        }

        public bool Modified
        {
            get { return modified; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public int ID
        {
            get { return 999; }
        }

        public void ClearState()
        {
            // possibly TODO 
            //if D_tan of initial state is wanted copy calculations from update material for 
            // d_prev_step=0 kai Delta [i] = 0 gia i=0,1 kai 2
            //but don't use it as elastic in other cases
            // maybe in iterative procedure (example provider.Reset ?)
        }

        public void ClearStresses() => throw new NotImplementedException();

        public void ClearTractions()
        {

        }

        public double YoungModulus
        {
            get { return youngModulus; }
            set { throw new InvalidOperationException(); }
        }

        public double PoissonRatio
        {
            get { return poissonRatio; }
            set { throw new InvalidOperationException(); }
        }

        private double[] coordinates;

        public double[] Coordinates
        {

            get { return coordinates; }
            set { throw new InvalidOperationException(); }
        }
    }
}
