using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
//using ISAAR.MSolve.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Triangulation;

namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// cohesion only friction type response in shear separation mode (with kinematic hardening like behaviour). Linear elastic behaviour in normal mode
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class BondSlipCohMatUniaxial : ICohesiveZoneMaterial3D // TODOGerasimos
    {
        private bool modified; // opws sto MohrCoulomb gia to modified

        public double k_elastic { get; set; } // opws sto elastic 3d 
        public double k_elastic2 { get; set; }
        public double k_elastic_normal { get; set; }
        public double t_max { get; set; }
        public double[] s_0 { get; set; }
        public double[] a_0 { get; set; }
        private double[] alpha { get; set; }
        public double tol { get; set; }
        private double[] eLastConverged;
        private double[] eCurrentUpdate;
        private double[,] ConstitutiveMatrix3D;
        private double[,] ConstitutiveMatrix3Dprevious;
        private BondSlipCohMat slipMaterial;
        private double[] stress3D;

        public BondSlipCohMatUniaxial(double k_elastic, double k_elastic2, double k_elastic_normal, double t_max, double[] s_0, double[] a_0, double tol)
        {
            this.k_elastic = k_elastic;
            this.k_elastic2 = k_elastic2;
            this.k_elastic_normal = k_elastic_normal;
            this.t_max = t_max;
            this.s_0 = s_0; //length = 2
            this.a_0 = a_0; //length = 2
            this.tol = tol;
            this.slipMaterial = new BondSlipCohMat(k_elastic, k_elastic2, k_elastic_normal, t_max, s_0, a_0, tol);
            this.InitializeMatrices();
        }

        public BondSlipCohMatUniaxial(double T_o_1, double D_o_1, double k_elastic2_ratio, double T_o_3, double D_o_3, double[] s_0, double[] a_0, double tol)
        {
            this.k_elastic = T_o_1 / D_o_1;
            this.k_elastic2 = k_elastic2_ratio * k_elastic;
            this.k_elastic_normal = T_o_3 / D_o_3;
            this.t_max = T_o_1; // Prosoxh exei lifthei idio koino orio diarrohs fy sunolika gia th sunistamenh ths paramorfwshs anexarthtws dieftunshs
            this.s_0 = s_0; //length = 2
            this.a_0 = a_0; //length = 2
            this.tol = tol;
            this.slipMaterial = new BondSlipCohMat(T_o_1, D_o_1, k_elastic2_ratio, T_o_3, D_o_3, s_0, a_0, tol);
            this.InitializeMatrices();
        }

        //ICohesiveZoneMaterial3D ICohesiveZoneMaterial3D.Clone()
        //{
        //    return this.Clone();
        //}

        object ICloneable.Clone() => Clone();
        ICohesiveZoneMaterial3D ICohesiveZoneMaterial3D.Clone() => Clone();

        public BondSlipCohMatUniaxial Clone()
        {
            return new BondSlipCohMatUniaxial(k_elastic, k_elastic2, k_elastic_normal, t_max, s_0, a_0, tol);
        }

        private double c1;
        //private double[] sigma;

        private bool matrices_not_initialized = true;

        public void InitializeMatrices()
        {
            eCurrentUpdate = new double[2];
            eLastConverged = new double[2];// TODO: na ginetai update sto save state ennoeitai mazi me ta s_0 klp. Mporei na xrhsimopooithei h grammh apo thn arxh tou update material
            stress3D = new double[3];
            ConstitutiveMatrix3D = new double[3, 3];

            ConstitutiveMatrix3D[0, 0] = k_elastic; ConstitutiveMatrix3D[1, 1] = k_elastic_normal; ConstitutiveMatrix3D[2, 2] = k_elastic_normal;
            ConstitutiveMatrix3Dprevious = new double[3, 3];
            ConstitutiveMatrix3Dprevious[0, 0] = k_elastic; ConstitutiveMatrix3Dprevious[1, 1] = k_elastic_normal; ConstitutiveMatrix3Dprevious[2, 2] = k_elastic_normal;

            matrices_not_initialized = false;
        }

        public void UpdateMaterial(double[] epsilon)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int j = 0; j < 3; j++)
                {
                    ConstitutiveMatrix3Dprevious[k, j] = ConstitutiveMatrix3D[k, j];
                }
            }
            slipMaterial.UpdateMaterial(new double[3] { epsilon[0], 0, 0 });
            ConstitutiveMatrix3D[0, 0] = slipMaterial.ConstitutiveMatrix[0, 0];
            ConstitutiveMatrix3D[1, 1] = k_elastic_normal;
            ConstitutiveMatrix3D[2, 2] = k_elastic_normal;
            stress3D = new double[3] { slipMaterial.Tractions[0], k_elastic_normal * epsilon[1], k_elastic_normal * epsilon[2] };

            this.modified = CheckIfConstitutiveMatrixChanged();
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    if (Math.Abs(ConstitutiveMatrix3Dprevious[i, j] - ConstitutiveMatrix3D[i, j]) > 1e-10)
                        return true;

            return false;
        }

        public double[] Tractions // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return stress3D; }
        }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {

                return Matrix.CreateFromArray(ConstitutiveMatrix3D);
            }
        }

        public void SaveState()
        {
            slipMaterial.SaveState();
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
            get { return 1000; }
        }

        public void ClearState() // pithanws TODO 
        {
            //ean thelei to D_tan ths arxikhs katastashs tha epistrepsoume const me De
            // alla oxi ia na to xrhsimopoihsei gia elastiko se alles periptwseis
            //opws
            // sthn epanalhptikh diadikasia (opws px provider.Reset pou sumvainei se polles epanalipseis?)
        }
        public void ClearTractions()
        {

        }

        private double youngModulus = 1;
        public double YoungModulus
        {
            get { throw new InvalidOperationException(); }
            set { throw new InvalidOperationException(); }
        }

        private double poissonRatio = 1;
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

        public void ClearStresses()
        {
            throw new NotImplementedException();
        }
    }
}
