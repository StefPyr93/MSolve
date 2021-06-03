using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using Accord.Math.Decompositions;
using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Stochastic
{
    public class MazarsConcreteMaterial_v2 : IIsotropicContinuumMaterial3D
    {
        private bool modified;

        public double youngModulus { get; set; }
        public double poissonRatio { get; set; }
        public double Strain_0 { get; set; }
        public double At { get; set; }
        public double Bt { get; set; }
        public double Ac { get; set; }
        public double Bc { get; set; }
        public double Veta { get; set; }

        object ICloneable.Clone() => Clone();
        //IContinuumMaterial3D IContinuumMaterial3D.Clone() => Clone();
        public MazarsConcreteMaterial_v2 Clone()
        {
            return new MazarsConcreteMaterial_v2()
            {
                youngModulus = this.youngModulus,
                poissonRatio = this.poissonRatio,
                Strain_0 = this.Strain_0,
                At = this.At,
                Bt = this.Bt,
                Ac = this.Ac,
                Bc = this.Bc,
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
        private double[,] strain_t;
        private double[,] strain_c;
        private double sum_princ_strains;
        private double[] princ_strain_t;
        private double[] princ_strain_c;
        private double at;
        private double ac;
        private double dt;
        private double dc;
        private double[,] D_tan_f;
        private double tol;

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
            strain_t = new double[3, 3];
            strain_c = new double[3, 3];
            princ_strain_t = new double[3];
            princ_strain_c = new double[3];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
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
            strain_eq = 0.0;
            var strain_mat = new double[3, 3] { { strain[0], 0.5 * strain[3], 0.5 * strain[5] }, { 0.5 * strain[3], strain[1], 0.5 * strain[4] }, { 0.5 * strain[5], 0.5 * strain[4], strain[2] } };
            var eigen_strain_mat = new EigenvalueDecomposition(strain_mat);
            var sum_principal_strains = 0.0;
            for (int i = 0; i < 3; i++)
            {
                if (eigen_strain_mat.DiagonalMatrix[i, i] > 0)
                {
                    strain_eq += Math.Pow(eigen_strain_mat.DiagonalMatrix[i, i],2);
                    sum_principal_strains += eigen_strain_mat.DiagonalMatrix[i, i];
                }
            }
            strain_eq = Math.Sqrt(strain_eq);
            if (strain_eq <= strain_eq_prev) //&& dmg==dmg_prev
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
            else
            {
                for (int i = 0; i < 6; i++)
                {
                    stress_eff[i] = 0;
                    for (int j = 0; j < 6; j++)
                    {
                        stress_eff[i] += this.constitutiveMatrix[i, j] * strain[j];
                    }
                }
                var stress_mat = new double[3, 3] { { stress_eff[0], stress_eff[3], stress_eff[5] },
                                 { stress_eff[3], stress_eff[1], stress_eff[4] }, { stress_eff[5], stress_eff[4], stress_eff[2] } };
                var eigen_stress = new EigenvalueDecomposition(stress_mat);
                var eigen_stress_pos = new double[3, 3];
                var eigen_stress_neg = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    if (eigen_stress.DiagonalMatrix[i, i] > 0)
                        eigen_stress_pos[i, i] = eigen_stress.DiagonalMatrix[i, i];
                    else
                        eigen_stress_neg[i, i] = eigen_stress.DiagonalMatrix[i, i];
                }
                var stress_t_temp1 = new double[3, 3];
                var stress_c_temp1 = new double[3, 3];
                var stress_t_temp2 = new double[3, 3];
                var stress_c_temp2 = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        strain_t[i, j] = 0;
                        strain_c[i, j] = 0;
                        for (int k = 0; k < 3; k++)
                        {

                            stress_t_temp1[i, j] += eigen_stress_pos[i, k] * eigen_stress.Eigenvectors[j, k];
                            stress_c_temp1[i, j] += eigen_stress_neg[i, k] * eigen_stress.Eigenvectors[j, k];
                        }
                    }
                }
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            stress_t_temp2[i, j] += eigen_stress.Eigenvectors[i, k] * stress_t_temp1[k, j];
                            stress_c_temp2[i, j] += eigen_stress.Eigenvectors[i, k] * stress_c_temp1[k, j];
                        }
                    }
                }
                var strain_t_temp = new double[6];
                var strain_c_temp = new double[6];
                for (int i = 0; i < 6; i++)
                {
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 0] * stress_t_temp2[0, 0];
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 1] * stress_t_temp2[1, 1];
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 2] * stress_t_temp2[2, 2];
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 3] * stress_t_temp2[0, 1];
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 4] * stress_t_temp2[1, 2];
                    strain_t_temp[i] += inverseConstitutiveMatrix[i, 5] * stress_t_temp2[0, 2];

                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 0] * stress_c_temp2[0, 0];
                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 1] * stress_c_temp2[1, 1];
                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 2] * stress_c_temp2[2, 2];
                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 3] * stress_c_temp2[0, 1];
                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 4] * stress_c_temp2[1, 2];
                    strain_c_temp[i] += inverseConstitutiveMatrix[i, 5] * stress_c_temp2[0, 2];
                }
                strain_t = new double[3, 3] { { strain_t_temp[0], strain_t_temp[3], strain_t_temp[5] }, { strain_t_temp[3], strain_t_temp[1], strain_t_temp[4] }, { strain_t_temp[5], strain_t_temp[4], strain_t_temp[2] } };
                strain_c = new double[3, 3] { { strain_c_temp[0], strain_c_temp[3], strain_c_temp[5] }, { strain_c_temp[3], strain_c_temp[1], strain_c_temp[4] }, { strain_c_temp[5], strain_c_temp[4], strain_c_temp[2] } };
                var eigen_strain_t = new EigenvalueDecomposition(strain_t);
                var eigen_strain_c = new EigenvalueDecomposition(strain_c);
                for (int i = 0; i < 3; i++)
                {
                    princ_strain_t[i] = eigen_strain_t.DiagonalMatrix[i, i];
                    princ_strain_c[i] = eigen_strain_c.DiagonalMatrix[i, i];
                }
                sum_princ_strains = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    if (princ_strain_t[i] > 0)
                    {
                        sum_princ_strains += princ_strain_t[i];
                    }
                    if (princ_strain_c[i] > 0)
                    {
                        sum_princ_strains += princ_strain_c[i];
                    }
                }
                at = 0;
                ac = 0;
                dt = 0;
                dc = 0;
                for (int i = 0; i < 3; i++)
                {
                    if (princ_strain_t[i] > 0)
                    {
                        at += Math.Pow(princ_strain_t[i] / sum_princ_strains, Veta);
                    }
                    if (princ_strain_c[i] > 0)
                    {
                        ac += Math.Pow(princ_strain_c[i] / sum_princ_strains, Veta);
                    }
                }
                //ac = 1;at = 0;
                dt = 1 - (1 - At) * Strain_0 / strain_eq - At / Math.Exp(Bt * (strain_eq - Strain_0));
                dc = 1 - (1 - Ac) * Strain_0 / strain_eq - Ac / Math.Exp(Bc * (strain_eq - Strain_0));
                dmg = at * dt + ac * dc;
                if (dmg < 0)
                    dmg = 0;
                else if (dmg > 1)
                    dmg = 1;
                var Dtan_coeff = at * ((1 - At) * Strain_0 / Math.Pow(strain_eq, 2) + At * Bt * Math.Exp(-Bt * (strain_eq - Strain_0))) +
                                  ac * ((1 - Ac) * Strain_0 / Math.Pow(strain_eq, 2) + Ac * Bc * Math.Exp(-Bc * (strain_eq - Strain_0))); //sum_principal_strains

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        if (strain[j] > 0)
                            D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j] - (Dtan_coeff * stress_eff[i]) * (strain[j] / strain_eq);
                        else
                            D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j];
                        D_tan_f[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j];
                        //D_tan[i, j] = (1 - dmg) * this.constitutiveMatrix[i, j]; //valid only for secant stiffness method
                    }
                }
                for (int i = 0; i < 6; i++)
                {
                    stress[i] = 0;
                    for (int j = 0; j < 6; j++)
                    {
                        stress[i] += D_tan_f[i, j] * strain[j];
                    }
                }
            }
            if(double.IsNaN(D_tan[0,0])) throw new ArithmeticException();
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
            double fE2 = fE1 * - PoissonRatio;
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
            if(strain_eq > strain_eq_prev)
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
