using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeperateIntegrationClassCheck
    {
        public static void RunExample() //(Model, double[]) RunExample()
        {
            IContinuumMaterial3DDefGrad material1 = new ElasticMaterial3DDefGrad() { PoissonRatio = 0.4, YoungModulus = 3.5 };
            double[,] consMatrix1 = material1.ConstitutiveMatrix.CopytoArray2D();

            int subdiscr1 = 4;// 6;
            int discr1 = 3;//4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = discr1*subdiscr1;// 23;
            int subdiscr1_shell = 6;//14;
            int discr1_shell = 1;
            int graphene_sheets_number = 3; //periektikothta 0.525% 

            double scale_factor = 2;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            

            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelopHSTAM(1, true, mpgp,
            //subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);
            //IContinuumMaterial3DDefGrad microstructure2 = new MicrostructureDefGrad3D(rveBuilder, true, 1);
            //double[,] consMatrix2 = microstructure2.ConstitutiveMatrix.CopytoArray2D();


        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            //// DUPLICATE ANNY CHANGES IN ONERVEEXAMPLEMPI
            //// kai sto integration Separatecodecheckingclass....GitTest.cs

            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            {
                rveMatrixParameters mp;
                mp = new rveMatrixParameters()
                {
                    E_disp = 3.5, //Gpa
                    ni_disp = 0.4, // stather Poisson
                    L01 = 95, //150, // diastaseis
                    L02 = 95, //150,
                    L03 = 95, //40,
                    hexa1 = discr1 * subdiscr1,// diakritopoihsh
                    hexa2 = discr1 * subdiscr1,
                    hexa3 = discr1 * subdiscr1,
                };

                grapheneSheetParameters gp;
                gp = new grapheneSheetParameters()
                {
                    // parametroi shell
                    E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                    ni_shell = 0.0607, // stathera poisson
                    elem1 = discr1_shell * subdiscr1_shell,
                    elem2 = discr1_shell * subdiscr1_shell,
                    L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                    L2 = 50,// nm
                    L3 = 112.5096153846, // nm
                    a1_shell = 0, // nm
                    tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                    //parametroi cohesive epifaneias
                    T_o_3 = 0.20, //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                    D_o_3 = 0.25, //0.5, // nm
                    D_f_3 = 4, // nm
                    T_o_1 = 0.20, //0.05,// Gpa
                    D_o_1 = 0.25, //0.5, // nm
                    D_f_1 = 4, // nm
                    n_curve = 1.4
                };

                Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
                return gpmp;
            }

            if (CnstValues.parameterSet == ParameterSet.stiffLargerRve)
            {

                rveMatrixParameters mp;
                mp = new rveMatrixParameters()
                {
                    E_disp = 3.5, //Gpa
                    ni_disp = 0.4, // stather Poisson
                    L01 = 120, //95, //150, // diastaseis
                    L02 = 120, //95, //150,
                    L03 = 120, //95, //40,
                    hexa1 = discr1 * subdiscr1,// diakritopoihsh
                    hexa2 = discr1 * subdiscr1,
                    hexa3 = discr1 * subdiscr1,
                };

                grapheneSheetParameters gp;
                gp = new grapheneSheetParameters()
                {
                    // parametroi shell
                    E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                    ni_shell = 0.0607, // stathera poisson
                    elem1 = discr1_shell * subdiscr1_shell,
                    elem2 = discr1_shell * subdiscr1_shell,
                    L1 = 38, //50,// nm  // DIORTHOSI 2 graphene sheets
                    L2 = 38, //.50,// nm
                    L3 = 112.5096153846, // nm
                    a1_shell = 0, // nm
                    tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                    //parametroi cohesive epifaneias
                    T_o_3 = 0.20, //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                    D_o_3 = 0.25, //0.5, // nm
                    D_f_3 = 4, // nm
                    T_o_1 = 0.20, //0.05,// Gpa
                    D_o_1 = 0.25, //0.5, // nm
                    D_f_1 = 4, // nm
                    n_curve = 1.4
                };

                Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
                return gpmp;
            }
            else { throw new NotImplementedException(); }
        }
    }
}
