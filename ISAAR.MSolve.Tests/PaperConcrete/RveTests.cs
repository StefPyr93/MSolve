﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.FEM.Materials;

using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MSAnalysis.RveTemplatesPaper;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
//using ISAAR.MSolve.MSAnalysis.RveTemplatesPaper;
//using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Tests.FEMpartB;
using MGroup.Stochastic;
//using JetBrains.dotMemoryUnit;
using MathNet.Numerics.LinearAlgebra;
using Troschuetz.Random;
using Xunit;
using MatlabWriter = MathNet.Numerics.Data.Matlab.MatlabWriter;

//[assembly: SuppressXUnitOutputException]

namespace ISAAR.MSolve.Tests
{
    public class RveTests
    {
        [Fact]
        public static void Check3DscaleTransitionsAndMicrostructure()
        {
            //MATERIAL LAW TO TEST WITH
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[6] { 0.01, 0, 0, 0, 0, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck1 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };
            //material1.SaveState();
            GLVec = new double[6] { 0, 0, 0, 0, 0.02, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck2 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };


            //TWO PHASE RVE BUILDER
            var outterMaterial = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            var innerMaterial = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };


            //WITH GMSH GEOMETRY DATA
            var homogeneousRveBuilder1 =
                new GmshCompositeRveBuilder(outterMaterial, innerMaterial, 2, 2, 2, "..\\..\\..\\RveTemplates\\Input\\Continuum\\1ball.msh");


            IContinuumMaterial3D microstructure3 = new Microstructure3D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[6] { 0.010, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[6] { 0, 0, 0, 0, 0.020, 0 });
            double[] stressesCheck4 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };

            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[6] { 0.030, 0, 0, 0, 0, 0 });
            double[] stressesCheck5 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            //COMPARISON
            Assert.True(AreDisplacementsSame(stressesCheck1, stressesCheck3));
            Assert.True(AreDisplacementsSame(stressesCheck2, stressesCheck4));
            Assert.True(AreDisplacementsSame(new double[6] { 3 * stressesCheck1[0], 3 * stressesCheck1[1], 3 * stressesCheck1[2], 3 * stressesCheck1[3], 3 * stressesCheck1[4], 3 * stressesCheck1[5] },
                                                                            stressesCheck5));
            Assert.True(AreDisplacementsSame(consCheck1, material1.ConstitutiveMatrix));
            Assert.True(AreDisplacementsSame(Matrix1.CopyToArray2D(), material1.ConstitutiveMatrix));
        }

        public static bool AreDisplacementsSame(double[,] expectedValues,
            IMatrixView computedValues)
        {
            var comparer = new ValueComparer(1E-14);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public static bool AreDisplacementsSame(double[] expectedValues,
            double[] computedValues, double tol = 1E-14)
        {
            var comparer = new ValueComparer(tol);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
                {
                    return false;
                }
            }
            return true;
        }

        [Fact]
        public void TestElasticAndMultiscaleMatricesTet()
        {
            var material3 = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210,
                PoissonRatio = 0.0,
                TangentVectorV1 = new double[3] { 1, 0, 0 },
                TangentVectorV2 = new double[3] { 0, 1, 0 }
            };

            //var trandom = new TRandom();
            //var randomInnerE = trandom.Normal(3.4e9, 0.2e9);
            var outterMaterial = new ElasticMaterial3DTotalStrain(0, 4.3210);

            var innerMaterial = new ElasticMaterial3DTotalStrain(0, 34);


            var homogeneousRveBuilder1 =
                new CompositeMaterialModeluilderTet2(outterMaterial, innerMaterial, 100, 100, 100);

            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] { 1, 0, 0 },
                TangentVectorV2 = new double[3] { 0, 1, 0 }
            };
            material4.UpdateMaterial(new double[] { 0, 0, 0 });
        }

        [Fact]
        public void CompositeRveWithElastic()
        {
            var material3 = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210,
                PoissonRatio = 0.0,
                TangentVectorV1 = new double[3] { 1, 0, 0 },
                TangentVectorV2 = new double[3] { 0, 1, 0 }
            };

            //var trandom = new TRandom();
            //var randomInnerE = trandom.Normal(3.4e9, 0.2e9);
            var outterMaterial = new ElasticMaterial3DTotalStrain(0, 4.3210);

            var innerMaterial = new ElasticMaterial3DTotalStrain(0, 4.3210);


            var homogeneousRveBuilder1 =
                new GmshCompositeRveBuilder(outterMaterial, innerMaterial, 2, 2, 2, "..\\..\\..\\RveTemplates\\Input\\Continuum\\t16Solid_physical_entities_no_volume_tag_change_More_inclusions.msh");

            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), true, 1)
            {
                TangentVectorV1 = new double[3] { 1, 0, 0 },
                TangentVectorV2 = new double[3] { 0, 1, 0 }
            };
            //material4.UpdateMaterial(new double[] { 0, 0, 0 });

            var cons = material4.ConstitutiveMatrix;

            Assert.Equal(material3.ConstitutiveMatrix[0, 0], cons[0, 0], 7);
            Assert.Equal(material3.ConstitutiveMatrix[1, 1], cons[1, 1], 7);
            Assert.Equal(material3.ConstitutiveMatrix[2, 2], cons[2, 2], 7);
        }

        [Fact]
        public static void CheckRVEwithCNTs()
        {
            //MATERIAL LAW TO TEST WITH
            double E_disp = 4; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[6] { 0.01, 0, 0, 0, 0, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck1 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };
            //material1.SaveState();
            GLVec = new double[6] { 0, 0, 0, 0, 0.02, 0 };
            material1.UpdateMaterial(GLVec);
            double[] stressesCheck2 = new double[6] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2], material1.Stresses[3], material1.Stresses[4], material1.Stresses[5] };


            //RVE BUILDER
            var matrixMaterial = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            int numberOfCnts = 1;

            var homogeneousRveBuilder1 =
                new CntReinforcedElasticNanocomposite(numberOfCnts);

            IContinuumMaterial3D microstructure3 = new Microstructure3D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            for (int i = 1; i <= 10; i++)
            {
                microstructure3.UpdateMaterial(new double[6] { 0.001 * i, 0, 0, 0, 0, 0 });
                double[] stressesCheck3 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
                var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
                microstructure3.SaveState();
            }

            microstructure3.UpdateMaterial(new double[6] { 0, 0, 0, 0, 0.020, 0 });
            double[] stressesCheck4 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };

            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[6] { 0.050, 0, 0, 0, 0, 0 });
            double[] stressesCheck5 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
            var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            //COMPARISON
            //Assert.True(AreDisplacementsSame(stressesCheck1, stressesCheck3));
            Assert.True(AreDisplacementsSame(stressesCheck2, stressesCheck4));
            Assert.True(AreDisplacementsSame(new double[6] { 3 * stressesCheck1[0], 3 * stressesCheck1[1], 3 * stressesCheck1[2], 3 * stressesCheck1[3], 3 * stressesCheck1[4], 3 * stressesCheck1[5] },
                                                                            stressesCheck5));
            Assert.True(AreDisplacementsSame(consCheck1, material1.ConstitutiveMatrix));
            Assert.True(AreDisplacementsSame(Matrix2.CopyToArray2D(), material1.ConstitutiveMatrix));
        }

        [Fact]
        public static void GenerateRVEwithCNTsSolutions()
        {
            LinearAlgebra.LibrarySettings.LinearAlgebraProviders = LinearAlgebra.LinearAlgebraProviderChoice.MKL;

            string pathName = @"C:\Users\stefp\OneDrive\Desktop\Solution Soimiri\matlabGeneratedCNTs\RVE Solutions";

            string InputFileName = "Input data.txt";
            string InputExtension = Path.GetExtension(InputFileName);
            string InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(InputFileName));
            string inputFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);

            string OutputFileName = "Output data.txt";
            string OutputExtension = Path.GetExtension(OutputFileName);
            string OutputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(OutputFileName));
            string outputFile = string.Format("{0}{1}", OutputfileNameOnly, OutputExtension);

            bool append = false;

            int numberOfCnts = 780;
            int solutions = 1;
            int increments_per_solution = 10;
            double[][] Input = new double[solutions * increments_per_solution][];
            double[][] Output = new double[solutions * increments_per_solution][];

            var homogeneousRveBuilder1 =
                new CntReinforcedElasticNanocomposite(numberOfCnts); //{ K_el = 20, K_pl = 2, T_max = 0.2, };
            homogeneousRveBuilder1.readFromText = true;

            IContinuumMaterial3D microstructure3 = new Microstructure3D(homogeneousRveBuilder1,
                model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

            for (int num_solution = 0; num_solution < solutions; num_solution++)
            {
                var MacroStrain = new double[6] {0.04, 0.02, 0.02, 0.04, 0.02, 0.02 };
                //var trandom = new TRandom();
                //for (int ii = 0; ii < 6; ii++) { MacroStrain[ii] = trandom.ContinuousUniform(-0.1, 0.1); }

                //homogeneousRveBuilder1.K_el = trandom.ContinuousUniform(0.1, 20);
                //homogeneousRveBuilder1.K_pl = trandom.ContinuousUniform(0.01, 2);
                //homogeneousRveBuilder1.T_max = trandom.ContinuousUniform(0.001, 0.2);

                //homogeneousRveBuilder1.UpdateCohesiveMaterial();

                //microstructure3 = new Microstructure3D(homogeneousRveBuilder1,
                //model => (new SuiteSparseSolver.Builder()).BuildSolver(model), false, 1);

                for (int i = 0; i < increments_per_solution; i++)
                {
                    var IncrMacroStrain = new double[6];
                    for (int ii = 0; ii < 6; ii++) { IncrMacroStrain[ii] = MacroStrain[ii] * (i + 1) / increments_per_solution; }
                    microstructure3.UpdateMaterial(new double[6] { IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] });
                    //Debug.WriteLine($"Strain {IncrMacroStrain[0]},{IncrMacroStrain[1]},{IncrMacroStrain[2]},{IncrMacroStrain[3]},{IncrMacroStrain[4]},{IncrMacroStrain[5]}");

                    double[] IncrMacroStress = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };

                    microstructure3.SaveState();
                    Input[num_solution * increments_per_solution + i] = new double[9] { homogeneousRveBuilder1.K_el, homogeneousRveBuilder1.K_pl, homogeneousRveBuilder1.T_max,
                                                                                        IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] };
                    Output[num_solution * increments_per_solution + i] = new double[6] { IncrMacroStress[0], IncrMacroStress[1], IncrMacroStress[2], IncrMacroStress[3], IncrMacroStress[4], IncrMacroStress[5] };

                    using (var writer = new StreamWriter(inputFile, append)) // append mode to continue from previous increment
                    {
                        writer.WriteLine($"{Input[num_solution * increments_per_solution + i][0]}, {Input[num_solution * increments_per_solution + i][1]}, {Input[num_solution * increments_per_solution + i][2]}, " +
                            $"{Input[num_solution * increments_per_solution + i][3]}, {Input[num_solution * increments_per_solution + i][4]}, {Input[num_solution * increments_per_solution + i][5]}, {Input[num_solution * increments_per_solution + i][6]}, " +
                            $"{Input[num_solution * increments_per_solution + i][7]}, {Input[num_solution * increments_per_solution + i][8]}");
                    }

                    using (var writer = new StreamWriter(outputFile, append)) // append mode to continue from previous increment
                    {
                        writer.WriteLine($"{Output[num_solution * increments_per_solution + i][0]}, {Output[num_solution * increments_per_solution + i][1]}, {Output[num_solution * increments_per_solution + i][2]}, " +
                            $"{Output[num_solution * increments_per_solution + i][3]}, {Output[num_solution * increments_per_solution + i][4]}, {Output[num_solution * increments_per_solution + i][5]}");
                    }
                    append = true;
                }
            }
        }

        [Fact]
        public static void TestNeuralNetwork()
        {
            LinearAlgebra.LibrarySettings.LinearAlgebraProviders = LinearAlgebra.LinearAlgebraProviderChoice.MKL;

            var NNtest = new NeuralNetwork();
            NNtest.ExtractNeuralNetworkParametersFromMatlab();

            var input = new double[6] { 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 };
            NNtest.constParameters = new double[3] { 10, 1, 0.1 };
            NNtest.InitializePreprocessing();
            for (int i = 0; i < 1131144; i++)
            {
                var output = NNtest.CalculateNeuralNetworkOutput(input);
                //var jacobian = NNtest.CalculateNeuralNetworkJacobian(input);
            }
        }
    }
}
