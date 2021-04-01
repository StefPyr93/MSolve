using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
//using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Loading.NodalLoads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MSAnalysis.RveTemplatesPaper;
using ISAAR.MSolve.MultiscaleAnalysis;
//using ISAAR.MSolve.MSAnalysis.RveTemplatesPaper;
//using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
//using JetBrains.dotMemoryUnit;
using MathNet.Numerics.LinearAlgebra;
using Troschuetz.Random;
using Xunit;
using MatlabWriter = MathNet.Numerics.Data.Matlab.MatlabWriter;

//[assembly: SuppressXUnitOutputException]

namespace ISAAR.MSolve.IGA.Tests
{
    public class PaperExamples
    {
        //[Fact]
        //public static void ScordelisLoShell()
        //{
        //    var filename = "ScordelisLoShell";
        //    var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt")
        //        .ToString(CultureInfo.InvariantCulture);
        //    var numberOfRealizations = 1000;

        //    var trandom1 = new TRandom();
        //    var trandom2 = new TRandom();
        //    var trandom3 = new TRandom();
        //    Stopwatch watch= new Stopwatch();
        //    watch.Start();
        //    for (int it = 0; it < 20; it++)
        //    {
        //        var youngModulusSolutionPairs = new double[numberOfRealizations, 1];
        //    Parallel.For(0, numberOfRealizations, realization =>
        //    {
        //        #region mat1
        //        var randomOutterE1 = trandom1.Lognormal(21.9056, 0.2880);

        //        var outterMaterial1 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = randomOutterE1,
        //            PoissonRatio = 0.0
        //        };
        //        var innerMaterial1 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = 1,
        //            PoissonRatio = 0.0
        //        };

        //        var homogeneousRveBuilder1 =
        //            new CompositeMaterialModeluilderTet2(outterMaterial1, innerMaterial1, 100, 100, 100);

        //        var material1 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder1,
        //            constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));
        //        #endregion
        //        #region mat2
        //        var randomOutterE2 = trandom2.Lognormal(21.9056, 0.2880);

        //        var outterMaterial2 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = randomOutterE2,
        //            PoissonRatio = 0.0
        //        };
        //        var innerMaterial2 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = 1,
        //            PoissonRatio = 0.0
        //        };

        //        var homogeneousRveBuilder2 =
        //            new CompositeMaterialModeluilderTet2(outterMaterial2, innerMaterial2, 100, 100, 100);

        //        var material2 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder2,
        //            constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));
        //        #endregion
        //        #region mat3
        //        var randomOutterE3 = trandom3.Lognormal(21.9056, 0.2880);

        //        var outterMaterial3 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = randomOutterE3,
        //            PoissonRatio = 0.0
        //        };
        //        var innerMaterial3 = new ElasticMaterial3DtotalStrain()
        //        {
        //            YoungModulus = 1,
        //            PoissonRatio = 0.0
        //        };

        //        var homogeneousRveBuilder3 =
        //            new CompositeMaterialModeluilderTet2(outterMaterial3, innerMaterial3, 100, 100, 100);

        //        var material3 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder3,
        //            constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));
        //        #endregion

        //        var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filepath, materials:new List<IShellMaterial>()
        //        {
        //            material1,
        //            material2, 
        //            material3
        //        });
        //        var model = modelReader.GenerateModelFromFile();

        //        model.SurfaceLoads.Add(new SurfaceDistributedLoad(-90, StructuralDof.TranslationY));

        //        // Rigid diaphragm for AB
        //        for (var i = 0; i < 19; i++)
        //        {
        //            model.ControlPointsDictionary[i * 19].Constraints
        //                .Add(new Constraint() {DOF = StructuralDof.TranslationX});
        //            model.ControlPointsDictionary[i * 19].Constraints
        //                .Add(new Constraint() {DOF = StructuralDof.TranslationY});
        //        }

        //        // Symmetry for CD
        //        for (var i = 0; i < 19; i++)
        //        {
        //            model.ControlPointsDictionary[i * 19 + 18].Constraints
        //                .Add(new Constraint() {DOF = StructuralDof.TranslationZ});

        //            model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
        //                new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationX),
        //                new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationX)));
        //            model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
        //                new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationY),
        //                new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationY)));
        //        }

        //        // Symmetry for AD
        //        for (var j = 0; j < 19; j++)
        //        {
        //            model.ControlPointsDictionary[j].Constraints
        //                .Add(new Constraint() {DOF = StructuralDof.TranslationX});
        //            model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
        //                new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationY),
        //                new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationY)));
        //            model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
        //                new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationZ),
        //                new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationZ)));
        //        }

        //        // Solvers
        //        var solverBuilder = new SkylineSolver.Builder();
        //        ISolver solver = solverBuilder.BuildSolver(model);

        //        // Structural problem provider
        //        var provider = new ProblemStructural(model, solver);

        //        // Linear static analysis
        //        var childAnalyzer = new LinearAnalyzer(model, solver, provider);
        //        var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

        //        // Run the analysis
        //        parentAnalyzer.Initialize();
        //        parentAnalyzer.Solve();


        //        var cp = model.ControlPointsDictionary.Values.Last();
        //        var dofA = model.GlobalDofOrdering.GlobalFreeDofs[cp, StructuralDof.TranslationY];

        //        var solution = solver.LinearSystems[0].Solution[dofA];

        //        youngModulusSolutionPairs[realization, 0] = solution;
        //    });

        //    var writer = new Array2DWriter();
        //    writer.WriteToFile(youngModulusSolutionPairs,
        //        Path.Combine(Directory.GetCurrentDirectory(), $"ScordelisLoMultiscaleResultsThicknessLayers_{it}.txt"));
        //    }
        //    watch.Stop();
        //    TimeSpan ts = watch.Elapsed;
        //    string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
        //        ts.Hours, ts.Minutes, ts.Seconds,
        //        ts.Milliseconds / 10);
        //    Console.WriteLine("RunTime " + elapsedTime);
        //}



        ////[Fact]
        //public void SimpleHoodBenchmark()
        //{
        //    var filename = "attempt2";
        //    var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.iga");

        //    for (int d = 0; d < 1; d++)
        //    {
        //        var numberOfRealizations = 1;
        //        var youngModulusSolutionPairs = new double[numberOfRealizations, 2];

        //        for (int realization = 0; realization < numberOfRealizations; realization++)
        //        {
        //            var randomCnts = 800;
        //            youngModulusSolutionPairs[realization, 0] = randomCnts;
        //            var outterMaterial = new ElasticMaterial3DTotalStrain(0.4, 4);
                    
        //            var homogeneousRveBuilder1 = new CntReinforcedElasticNanocomposite(outterMaterial, randomCnts);

        //            var material4 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder1,
        //                constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));

        //            var model = new Model();
        //            var modelReader = new IgaFileReader(model, filepath);
        //            modelReader.CreateTSplineShellsModelFromFile(IgaFileReader.TSplineShellType.Thickness, material4);

        //            for (int i = 0; i < 100; i++)
        //            {
        //                var id = model.ControlPoints.ToList()[i].ID;
        //                model.ControlPointsDictionary[id].Constraints
        //                    .Add(new Constraint() {DOF = StructuralDof.TranslationX});
        //                model.ControlPointsDictionary[id].Constraints
        //                    .Add(new Constraint() {DOF = StructuralDof.TranslationY});
        //                model.ControlPointsDictionary[id].Constraints
        //                    .Add(new Constraint() {DOF = StructuralDof.TranslationZ});
        //            }

        //            for (int i = model.ControlPoints.Count() - 100; i < model.ControlPoints.Count(); i++)
        //            {
        //                var id = model.ControlPoints.ToList()[i].ID;
        //                model.Loads.Add(new NodalLoad(model.ControlPointsDictionary[id],StructuralDof.TranslationZ,100/10e4));
        //            }

        //            var solverBuilder = new SuiteSparseSolver.Builder();
        //            ISolver solver = solverBuilder.BuildSolver(model);

        //            // Structural problem provider
        //            var provider = new ProblemStructural(model, solver);

        //            // Linear static analysis
        //            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
        //            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

        //            // Run the analysis
        //            parentAnalyzer.Initialize();

        //            parentAnalyzer.Solve();
                    
        //            youngModulusSolutionPairs[realization, 1] = solver.LinearSystems[0].Solution[46448];
        //        }

        //        var writer = new Array2DWriter();
        //        writer.WriteToFile(youngModulusSolutionPairs,
        //            Path.Combine(Directory.GetCurrentDirectory(), $"BumperMultiscaleResults_{d}"));
        //    }

        //}

        ////[Fact]
        //public void TestElasticAndMultiscaleMatricesCnt()
        //{
        //    var numberOfRealizations = 20000;
        //    var youngModulusSolutionPairs = new double[numberOfRealizations, 2];
        //    for (int i = 0; i < numberOfRealizations; i++)
        //    {
        //        var trandom = new TRandom();
        //        var randomCnts = trandom.Normal(200, 20);

        //        var material3 = new ShellElasticMaterial2Dtransformationb()
        //        {
        //            YoungModulus = 4,
        //            PoissonRatio = 0.4,
        //            TangentVectorV1 = new double[3] {1, 0, 0},
        //            TangentVectorV2 = new double[3] {0, 1, 0}
        //        };

        //        var outterMaterial = new ElasticMaterial3DTotalStrain(0.4, 4);
                
        //        var homogeneousRveBuilder1 = new CntReinforcedElasticNanocomposite(outterMaterial, (int)randomCnts);

        //        var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
        //            model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
        //        {
        //            TangentVectorV1 = new double[3] {1, 0, 0},
        //            TangentVectorV2 = new double[3] {0, 1, 0}
        //        };
        //        material4.UpdateMaterial(new double[]{0,0,0});


        //        youngModulusSolutionPairs[i, 0] = (int)randomCnts;
        //        youngModulusSolutionPairs[i, 0] = material4.ConstitutiveMatrix[0,0];
        //    }
        //    var writer = new Array2DWriter();
        //    writer.WriteToFile(youngModulusSolutionPairs,
        //        Path.Combine(Directory.GetCurrentDirectory(), "materialCnt200"));
        //}


        //[Fact]
        public void TestElasticAndMultiscaleMatricesTet()
        {
            var material3 = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210,
                PoissonRatio = 0.0,
                TangentVectorV1 = new double[3] {1, 0, 0},
                TangentVectorV2 = new double[3] {0, 1, 0}
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
                TangentVectorV1 = new double[3] {1, 0, 0},
                TangentVectorV2 = new double[3] {0, 1, 0}
            };
            material4.UpdateMaterial(new double[]{0,0,0});
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


        ////[Fact]
        //public void BenchmarkEmbedding()
        //{
        //    for (int i = 0; i < 10; i++)
        //    {
        //        var outterMaterial = new ElasticMaterial3DTotalStrain(0.4, 4);

        //        var homogeneousRveBuilder1 = new CntReinforcedElasticNanocomposite(outterMaterial, 300);

        //        var material4 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder1,
        //            constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));
        //    }

        //}
    }
}