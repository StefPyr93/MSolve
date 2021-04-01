using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MSAnalysis.RveTemplates.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Solvers.Direct;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using Xunit;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class MemoryProfilerRveExample
    {
        public static (double[], double[])/* double[,], IVector, IVector)*/ CheckExample46InputInCodeRAMconsumption() //palio "Check_Graphene_rve_Obje_Integration()"
        {
            #region rve builder parameters and example choice
            CnstValues.exampleNo = 46; CnstValues.parameterSet = ParameterSet.stiffLargerRve;
            CnstValues.runOnlyHexaModel = false;
            CnstValues.PreventMATLABandTotalOutput();
            CnstValues.useInput_forRVE = true; //Panta prin thn getRveModelAndBoundaryNodes
            CnstValues.isInputInCode_forRVE = true;


            (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor) = SeparateCodeCheckingClass5b_bNEW_debugGit.GetGrRveExampleDiscrDataFromFile(new CnstValues());
            int discr3 = discr1 * subdiscr1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = SeperateIntegrationClassCheck.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            #endregion

            #region solve skyline Microstructures (with Git and GitSerial RveBuilders)                    
            //var rveBuilder3 = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            //subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);
            //var microstructure2Serial = new MicrostructureDefGrad3D(rveBuilder3,
            //    model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

            //double[,] consCheckSerial2 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheckSerial2[i1, i2] = microstructure2Serial.ConstitutiveMatrix[i1, i2]; } }
            //microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            //microstructure2Serial.SaveState();
            //microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            //Vector solutionSuiteSparse = (Vector)microstructure2Serial.uInitialFreeDOFDisplacementsPerSubdomain.ElementAt(0).Value.Copy();
            //double[] stressesSuitesparse = microstructure2Serial.Stresses;
            //double[,] constitutiveSuitesparse = microstructure2Serial.ConstitutiveMatrix.CopytoArray2D();
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveSuitesparse[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            #endregion

            #region solve microstructure with feti dp solver
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number, false);
            var microstructure3 = new MicrostructureDefGrad3DSerial(rveBuilder,
                rveBuilder.GetAppropriateSolverMpi, false, 1, true, true);


            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            microstructure3.SaveState();
            double[] sttressesFeti1 = new double[6];
            for (int i1 = 0; i1 < 6; i1++)
            {
                sttressesFeti1[i1] = microstructure3.Stresses[i1];
            }

            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesFeti = microstructure3.Stresses;
            double[,] constitutiveFeti = microstructure3.ConstitutiveMatrix.CopytoArray2D();
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveFeti[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.SaveState();
            //Vector SolutionFetiInSerialFormat = ReorderSolutionInSerialSkylineFormat(microstructure2Serial.model, solutionSuiteSparse, microstructure3.model, microstructure3.uInitialFreeDOFDisplacementsPerSubdomain);
            #endregion

            //var errorVector = (vector2 - globalUvectrInSerialFormat).Scale(1 / vector2.Norm2());
            //double errorNorm = (solutionSuiteSparse - SolutionFetiInSerialFormat).Norm2() / solutionSuiteSparse.Norm2();

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, stressesSuitesparse, 1e-8));
            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(solutionSuiteSparse.CopyToArray(), SolutionFetiInSerialFormat.CopyToArray(), 1e-4));

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(constitutiveSuitesparse, constitutiveFeti, 1e-7));


            //PrintUtilities.WriteToFileVector(stressesCheck3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck3.txt");
            //PrintUtilities.WriteToFile(consCheck1, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\consCheck1.txt");

            return (sttressesFeti1, stressesFeti);
        }

        public static (double[], double[])/*, double[,], IVector, IVector)*/  CheckExample46InputInCodeRAMconsumptionV2() //palio "Check_Graphene_rve_Obje_Integration()"
        {
            #region rve builder parameters and example choice
            CnstValues.exampleNo = 46; CnstValues.parameterSet = ParameterSet.stiffLargerRve;
            CnstValues.runOnlyHexaModel = false;
            CnstValues.PreventMATLABandTotalOutput();
            CnstValues.useInput_forRVE = true; //Panta prin thn getRveModelAndBoundaryNodes
            CnstValues.isInputInCode_forRVE = true;
            CnstValues.useV2FiniteElements = true;


            (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor) = SeparateCodeCheckingClass5b_bNEW_debugGit.GetGrRveExampleDiscrDataFromFile(new CnstValues());
            int discr3 = discr1 * subdiscr1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = SeperateIntegrationClassCheck.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            #endregion

            #region solve skyline Microstructures (with Git and GitSerial RveBuilders)                    
            //var rveBuilder3 = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            //subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);
            //var microstructure2Serial = new MicrostructureDefGrad3D(rveBuilder3,
            //    model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

            //double[,] consCheckSerial2 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheckSerial2[i1, i2] = microstructure2Serial.ConstitutiveMatrix[i1, i2]; } }
            //microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            //microstructure2Serial.SaveState();
            //microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            //Vector solutionSuiteSparse = (Vector)microstructure2Serial.uInitialFreeDOFDisplacementsPerSubdomain.ElementAt(0).Value.Copy();
            //double[] stressesSuitesparse = microstructure2Serial.Stresses;
            //double[,] constitutiveSuitesparse = microstructure2Serial.ConstitutiveMatrix.CopytoArray2D();
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveSuitesparse[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            #endregion

            #region solve microstructure with feti dp solver
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number, false);
            var microstructure3 = new MicrostructureDefGrad3DSerial(rveBuilder,
                rveBuilder.GetAppropriateSolverMpi, false, 1, true, true);


            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] sttressesFeti1 = new double[6];
            for (int i1 = 0; i1 < 6; i1++)
            {
                sttressesFeti1[i1] = microstructure3.Stresses[i1];
            }


            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesFeti = microstructure3.Stresses;
            double[,] constitutiveFeti = microstructure3.ConstitutiveMatrix.CopytoArray2D();
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveFeti[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.SaveState();
            //Vector SolutionFetiInSerialFormat = ReorderSolutionInSerialSkylineFormat(microstructure2Serial.model, solutionSuiteSparse, microstructure3.model, microstructure3.uInitialFreeDOFDisplacementsPerSubdomain);
            #endregion

            //var errorVector = (vector2 - globalUvectrInSerialFormat).Scale(1 / vector2.Norm2());
            //double errorNorm = (solutionSuiteSparse - SolutionFetiInSerialFormat).Norm2() / solutionSuiteSparse.Norm2();

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, stressesSuitesparse, 1e-8));
            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(solutionSuiteSparse.CopyToArray(), SolutionFetiInSerialFormat.CopyToArray(), 1e-4));

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(constitutiveSuitesparse, constitutiveFeti, 1e-7));


            //PrintUtilities.WriteToFileVector(stressesCheck3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck3.txt");
            //PrintUtilities.WriteToFile(consCheck1, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\consCheck1.txt");

            //return (stressesSuitesparse, stressesFeti, constitutiveSuitesparse, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);
            return (sttressesFeti1, stressesFeti);
        }
    }
}
