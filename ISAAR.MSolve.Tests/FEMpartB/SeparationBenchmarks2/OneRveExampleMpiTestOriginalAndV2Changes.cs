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

namespace ISAAR.MSolve.Tests.FEMpartB.SeparationBenchmarks2
{
    public class OneRveExampleMpiTestOriginalAndV2Changes // palio: "SeparateCodeCheckingClass4 "
    {
        [Fact]
        public static /*(double[], double[], double[,], IVector, IVector)*/ void Check_Graphene_rve_serial() //palio "Check_Graphene_rve_Obje_Integration()"
        {
            #region rve builder parameters and example choice
            CnstValues.exampleNo = 58;
            CnstValues.isInputInCode_forRVE = true;
            CnstValues.useInput_forRVE = true; //Panta prin thn getRveModelAndBoundaryNodes
            

            (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor) = GetGrRveExampleDiscrDataFromFile(new CnstValues());
            int discr3 = discr1 * subdiscr1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            #endregion

            #region solve skyline Microstructures (with Git and GitSerial RveBuilders)                    
            var rveBuilder3 = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            var microstructure2Serial = new MicrostructureDefGrad3D(rveBuilder3,
                model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

            //double[,] consCheckSerial2 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheckSerial2[i1, i2] = microstructure2Serial.ConstitutiveMatrix[i1, i2]; } }
            microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            microstructure2Serial.SaveState();
            microstructure2Serial.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            Vector solutionSuiteSparse = (Vector)microstructure2Serial.uInitialFreeDOFDisplacementsPerSubdomain.ElementAt(0).Value.Copy();
            double[] stressesSuitesparse = microstructure2Serial.Stresses;
            double[,] constitutiveSuitesparse = microstructure2Serial.ConstitutiveMatrix.CopytoArray2D();

            #region write serial solver output
            //var cnst = new CnstValues();
            //(new ISAAR.MSolve.LinearAlgebra.Output.Array1DWriter()).WriteToFile(solutionSuiteSparse.CopyToArray(), cnst.exampleOutputPathGen + @"\oneRveExampleTestSuitsparseSolution.txt");
            //(new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(microstructure2Serial.ConstitutiveMatrix.CopytoArray2D(), cnst.exampleOutputPathGen + @"\oneRveExampleTestSuitsparseConstitutive.txt");

            //(new ISAAR.MSolve.LinearAlgebra.Output.Array1DWriter()).WriteToFile(microstructure2Serial.Stresses, cnst.exampleOutputPathGen + @"\oneRveExampleTestSuitsparseStresses.txt");
            #endregion

            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveSuitesparse[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            #endregion

            #region solve microstructure with feti dp solver
            //CnstValues.useV2FiniteElements = true;
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            var microstructure3 = new MicrostructureDefGrad3DSerial(rveBuilder,
                rveBuilder.GetAppropriateSolverMpi, false, 1, true, true);
            

            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.03, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesFeti = microstructure3.Stresses;
            double[,] constitutiveFeti = microstructure3.ConstitutiveMatrix.CopytoArray2D();
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { constitutiveFeti[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.SaveState();
            Vector SolutionFetiInSerialFormat = ReorderSolutionInSerialSkylineFormat(microstructure2Serial.model, solutionSuiteSparse, microstructure3.model, microstructure3.uInitialFreeDOFDisplacementsPerSubdomain);
            #endregion

            //var errorVector = (vector2 - globalUvectrInSerialFormat).Scale(1 / vector2.Norm2());
            double errorNorm = (solutionSuiteSparse - SolutionFetiInSerialFormat).Norm2() / solutionSuiteSparse.Norm2();

            double pcg_tol = 1e-5;
            if (pcg_tol == 1e-8)
            {
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, stressesSuitesparse, 1e-8));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, new double[6] { 0.23960891245982407, 0.15427370248128061, 0.153671112977783, 0.0020627410983659151, -0.00071535481695788922, -0.0019255840674238048 }, 1e-14));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(solutionSuiteSparse.CopyToArray(), SolutionFetiInSerialFormat.CopyToArray(), 1e-4));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(constitutiveSuitesparse, constitutiveFeti, 1e-7));
            }
            else if (pcg_tol==1e-5)
            {
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, stressesSuitesparse, 1e-5));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesFeti, new double[6] {0.23960891248433011, 0.154273702183357, 0.15367111266616779, 0.0020627410561165154, -0.00071535523131603337, -0.0019255841128574466 }, 1e-14));
            }

            
            


            //PrintUtilities.WriteToFileVector(stressesCheck3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck3.txt");
            //PrintUtilities.WriteToFile(consCheck1, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\consCheck1.txt");
            //return (stressesSuitesparse, stressesFeti, constitutiveSuitesparse, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);
        }
        
        [Fact]
        public static /*(double[], double[], double[,], IVector, IVector)*/ void CheckExample46InputInCodeRAMconsumption() //palio "Check_Graphene_rve_Obje_Integration()"
        {
            #region rve builder parameters and example choice
            CnstValues.exampleNo = 46; CnstValues.parameterSet = ParameterSet.stiffLargerRve;
            CnstValues.runOnlyHexaModel = false;
            CnstValues.PreventMATLABandTotalOutput();
            CnstValues.useInput_forRVE = true; //Panta prin thn getRveModelAndBoundaryNodes
            CnstValues.isInputInCode_forRVE = true;


            (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor) = GetGrRveExampleDiscrDataFromFile(new CnstValues());
            int discr3 = discr1 * subdiscr1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            #endregion

            #region solve skyline Microstructures (with Git and GitSerial RveBuilders)                    
            var rveBuilder3 = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            var microstructure2Serial = new MicrostructureDefGrad3D(rveBuilder3,
                model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

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
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            var microstructure3 = new MicrostructureDefGrad3DSerial(rveBuilder,
                rveBuilder.GetAppropriateSolverMpi, false, 1, true, true);
            

            microstructure3.UpdateMaterial(new double[9] { /*1.10*/ 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
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
        }
        private static Vector ReorderSolutionInSerialSkylineFormat(Model skylineModel, Vector skylineSolution, Model fetimodel, Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain)
        {
            var freeNodesIds = skylineModel.GlobalDofOrdering.GlobalFreeDofs.GetRows().Select(x => x.ID);
            Vector globalU1_2 = Vector.CreateZero(skylineSolution.Length);
            foreach (int nodeID in freeNodesIds)
            {
                var fetiSubdomain = fetimodel.GetNode(nodeID).SubdomainsDictionary.ElementAt(0).Value;
                foreach (var dof in fetimodel.GlobalDofOrdering.GlobalFreeDofs.GetDataOfRow(fetimodel.GetNode(nodeID)).Keys)
                {
                    int skylineModelDofOrder = skylineModel.GlobalDofOrdering.GlobalFreeDofs[skylineModel.GetNode(nodeID), dof];
                    int fetisubdomainDofOrder = fetiSubdomain.FreeDofOrdering.FreeDofs[fetimodel.GetNode(nodeID), dof];

                    globalU1_2[skylineModelDofOrder] = uInitialFreeDOFDisplacementsPerSubdomain[fetiSubdomain.ID][fetisubdomainDofOrder];
                }
            }

            return globalU1_2;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            //// DUPLICATE ANNY CHANGES IN SEPARATEINTEGRATIONCLASSCHEC

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

        public static int[] ReadIntVector(string path)
        {

            var reader = new StreamReader(path);
            var lines = File.ReadLines(path).Count();
            int[] data = new int[lines];
            for (int i = 0; i < lines; ++i)
            {
                data[i] = Convert.ToInt32(reader.ReadLine());

            }
            reader.Close();
            reader.Dispose();
            return data;
        }

        private static (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor) GetGrRveExampleDiscrDataFromFile(CnstValues cnstValues)
        {
            //DUPLICATE CHANGES IN SAMPLE CONSOLE
            if (!CnstValues.isInputInCode_forRVE)
            {
                int[] discrData = ReadIntVector(cnstValues.exampleDiscrInputPathGen + @"\subdiscr1_discr1_ subdiscr1_shell_discr1_shell_graphene_sheets_number" + ".txt");
                double[] modelScaleFactor = MultiscaleAnalysis.SupportiveClasses.PrintUtilities.ReadVector(cnstValues.exampleDiscrInputPathGen + @"\modelScalingFactor" + ".txt");

                return (discrData[0], discrData[1], discrData[2], discrData[3], discrData[4], modelScaleFactor[0]);
            }
            else
            {
                (int[] discrData, double[] modelScaleFactor) = GeometryProviderForMpi.GetDiscrDataAndModelScaleFactor();

                //    = ReadIntVector(cnstValues.exampleDiscrInputPathGen + @"\subdiscr1_discr1_ subdiscr1_shell_discr1_shell_graphene_sheets_number" + ".txt");
                //double[] modelScaleFactor = MultiscaleAnalysis.SupportiveClasses.PrintUtilities.ReadVector(cnstValues.exampleDiscrInputPathGen + @"\modelScalingFactor" + ".txt");

                return (discrData[0], discrData[1], discrData[2], discrData[3], discrData[4], modelScaleFactor[0]);
            }

        }

        public static (int[], int[], int[]) Check_Graphene_rve_parallel() //palio "Check_Graphene_rve_Obje_Integration()"
        {
            //Origin h methodos Check_Graphene_rve_serial() tou parontos
            //Origin: SeparateCodeCheckingClass4.Check_Graphene_rve_Obje_Integration apo to branch: example/ms_development_nl_elements_merge
            //modifications: update kai tha xrhsimopoithei o GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostData 
            //o opoios exei kai antistoixo ddm: GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataDdm pou tha trexei akrivws apo katw
            //PROSOXH gia na elegxei kai h defterh iteration u_sunol_micro_2 prepei na valoume ston graphenebuilder Addgraphenesheet xwris to bondslip.

            //mporoun na ginoun delete:
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(GLVec);
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(GLVec);
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses;

            // den xreiazetai poia VectorExtensions.AssignTotalAffinityCount();
            var grapheneRveBuilder1 = new RveGrShMultipleSeparated(1);
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            // pros to paron
            var ModelAndNodes = grapheneRveBuilder1.GetModelAndBoundaryNodes();
            int[] hexaPrint = grapheneRveBuilder1.hexaPrint;
            int[] cohePrint = grapheneRveBuilder1.cohePrint;
            int[] shellPrint = grapheneRveBuilder1.shellPrint;
            return (hexaPrint, cohePrint, shellPrint);

            IContinuumMaterial3DDefGrad microstructure3 = new MicrostructureDefGrad3D(grapheneRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses;
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses;


        }

        #region transformation methods
        public static double[] Transform_DGtr_to_GLvec(double[,] DGtr)
        {
            double[,] GL = new double[3, 3];

            //
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0;
                    for (int p = 0; p < 3; p++)
                    {
                        GL[m, n] += DGtr[m, p] * DGtr[n, p];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                GL[m, m] += -1;
            }
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0.5 * GL[m, n];
                }
            }

            double[] GLvec = new double[6];
            //
            for (int m = 0; m < 3; m++)
            {
                GLvec[m] = GL[m, m];
            }
            GLvec[3] = 2 * GL[0, 1];
            GLvec[4] = 2 * GL[1, 2];
            GLvec[5] = 2 * GL[2, 0];

            return GLvec;
        }
        #endregion

        #region methodoi ths palaias SeparateCodeCheckingClass4
        //public static void Check05bStressIntegrationObjeIntegration()
        //{
        //    //Origin: SeparateCodeCheckingClass.Check05bStressIntegration
        //    //modifications: tha xrhsimopoithei h nea microstructure me obje kapoia subdomainCalculations

        //    double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
        //    ElasticMaterial3D material1 = new ElasticMaterial3D()
        //    { YoungModulus = E_disp, PoissonRatio = ni_disp, };
        //    double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        //    double[] GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
        //    material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
        //    //double[] stressesCheck1 = material1.Stresses;
        //    double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
        //        material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
        //    DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        //    GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
        //    material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
        //    material1.SaveState();
        //    double[] stressesCheck2 = material1.Stresses.Data;

        //    VectorExtensions.AssignTotalAffinityCount();
        //    IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27Hexa();
        //    //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

        //    IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObj(homogeneousRveBuilder1, false, 1);
        //    //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
        //    double[,] consCheck1 = new double[6, 6];
        //    for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

        //    microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
        //    double[] stressesCheck3 = microstructure3.Stresses.Data;
        //    microstructure3.SaveState();
        //    microstructure3.UpdateMaterial(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });
        //    double[] stressesCheck4 = microstructure3.Stresses.Data;
        //}

        //public static void Check05bStressIntegrationObje_Integration()
        //{
        //    //Origin: SeparateCodeCheckingClass.Check05bStressIntegration
        //    //modifications: tha xrhsimopoithei h nea microstructure me obje kapoia subdomainCalculations

        //    double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
        //    ElasticMaterial3D material1 = new ElasticMaterial3D()
        //    { YoungModulus = E_disp, PoissonRatio = ni_disp, };
        //    double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        //    double[] GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
        //    material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
        //    //double[] stressesCheck1 = material1.Stresses;
        //    double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
        //        material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
        //    DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        //    GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
        //    material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
        //    material1.SaveState();
        //    double[] stressesCheck2 = material1.Stresses.Data;

        //    VectorExtensions.AssignTotalAffinityCount();
        //    IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27Hexa();
        //    //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

        //    IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObj(homogeneousRveBuilder1, new SkylineSolver.Builder(), false, 1);
        //    //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
        //    double[,] consCheck1 = new double[6, 6];
        //    for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

        //    microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
        //    double[] stressesCheck3 = microstructure3.Stresses.Data;
        //    microstructure3.SaveState();
        //    microstructure3.UpdateMaterial(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });
        //    double[] stressesCheck4 = microstructure3.Stresses.Data;
        //}
        #endregion
    }
}
