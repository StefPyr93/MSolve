using System;
using System.Collections.Generic;
using System.Text;
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
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using System.IO;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;

namespace ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    public static class RveDataPrintMethods
    {
        public static void CheckOutputWriteFile()
        {
            string modelOutputPath_gen  = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2d_alte";
            int subdiscr1 = 6;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 23;
            int subdiscr1_shell = 14;
            int discr1_shell = 1;
            var mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            var mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
            var gp = mpgp.Item2;
            WriteModelDataOutput(modelOutputPath_gen, subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, mp, gp, new int[6] { 1, 2, 3, 4, 5, 6 }, 
                new double[2][] { new double[2] { 0.5, 0.79 }, new double[2] { 0.85, 0.97 } });
        }

        public static void WriteModelDataOutput(string modelOutputPath_gen, int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell,
            rveMatrixParameters mp,  grapheneSheetParameters gp, int[] kanonas_renumnering_2, double[][] o_xsunol_vectors, int[] sunol_nodes_numbering=null)
        {
            string modelDataOutputPath = modelOutputPath_gen + @"\montelo_graphene_generated.m";

            string[] matlabLines = new string[51] {"cycles="+1.ToString()+';',
            "%fileID = fopen('extracted_data2.txt');",
            "extr_data=textscan(fileID,'%f %f %f');",
            "fclose(fileID);",
            "a = size(extr_data{ 1,1},1);",
            "dg_prosth = zeros(a, 3);",
            "dg_prosth(:, 1) = extr_data{ 1,1}(:,:);% (:,:)->(:, 1)",
            "dg_prosth(:, 2) = extr_data{ 1,2};",
            "dg_prosth(:, 3) = extr_data{ 1,3};",
            "global subdiscr1 discr1 discr3 subdiscr1_shell discr1_shell",
            "subdiscr1="+subdiscr1.ToString()+";",
            "discr1 =" +discr1.ToString()+"; % DEN ALLAZEI",
            "% discr2 % den xrhsimopoieitai idia me discr1",
            "discr3 ="+ discr3.ToString()+";",
            "subdiscr1_shell ="+subdiscr1_shell.ToString()+";",
            "discr1_shell ="+discr1_shell.ToString()+";",
            " ",
            "%parametroi shell",
            "E="+gp.E_shell.ToString()+"; % GPa = 1000Mpa = 1000N / mm2",
            "ni ="+gp.ni_shell.ToString()+"; % stathera poisson",
            "elem2 ="+gp.elem1.ToString()+";",
            "elem1 ="+gp.elem2.ToString()+";%discr1_shell * subdiscr1_shell;",
            "L1 = "+gp.L1.ToString()+";% nm",
            "L2 = "+gp.L2.ToString()+";% nm",
            "L3 = "+gp.L3.ToString()+";% nm",
            "a1_shell = "+gp.a1_shell.ToString()+"% nm",
            "tk = 0.0125016478913782 * ones(elem1 * (3 * elem2 + 2) + 2 * elem2 + 1, 1); % 0.0125016478913782nm",
            " ",
            "% parametroi cohesive epifaneias",
            "% T_o_3,D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve",
            " T_o_3 ="+gp.T_o_3.ToString()+";% Gpa = 1000Mpa = 1000N / mm2",
            "D_o_3 = "+gp.D_o_3.ToString()+"; % nm",
            "D_f_3 = "+gp.D_f_3.ToString()+";% nm",
            " ",
            "T_o_1 ="+gp.T_o_1.ToString()+";% Gpa",
            "D_o_1 ="+gp.D_o_1.ToString()+";% nm",
            "D_f_1 ="+gp.D_f_1.ToString()+";% nm",
            " ",
            "n_curve ="+gp.n_curve.ToString()+";",
            " ",
            "% parametroi matrix rve",
            "E_disp ="+mp.E_disp.ToString()+";% Gpa",
            "ni_disp ="+mp.ni_disp.ToString()+";% stathera Poisson",
            " ",
            "L01 = "+mp.L01.ToString()+"; % nm",
            "L02 = "+mp.L02.ToString()+";",
            "L03 = "+mp.L03.ToString()+";",
            " ",
            "hexa1 ="+mp.hexa1.ToString()+";%discr1 * subdiscr1;",
            "hexa2 ="+mp.hexa1.ToString()+";discr1 * subdiscr1;",
            "hexa3 ="+mp.hexa1.ToString()+";discr3; ",
        };

             DdmCalculationsGeneral.WriteToFileStringArray(matlabLines, modelDataOutputPath);

            string kanonas_ren_2_path = modelOutputPath_gen + @"\kanonas_renum_2_generated.txt";
                DdmCalculationsGeneral.WriteToFileVector(kanonas_renumnering_2, kanonas_ren_2_path);

            string file_name = @"\o_xsunol_generated_gs_{0}.txt";

            DdmCalculationsGeneral.WriteToFileVectorsWithCounter(o_xsunol_vectors, modelOutputPath_gen, file_name);

            if(sunol_nodes_numbering!=null)
            {
                string renumbering_path = modelOutputPath_gen + @"\sunol_nodes_numbering.txt";
                DdmCalculationsGeneral.WriteToFileVector(sunol_nodes_numbering, renumbering_path);
            }
            
        }

        public static void WriteToFile(string[] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {

                writer.Write(array[i]);
                //writer.Write(' ');

                writer.WriteLine();
            }
            writer.Flush();



            writer.Dispose();
        }
    }
}
