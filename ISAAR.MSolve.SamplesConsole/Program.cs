using System;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed.Tests;
using ISAAR.MSolve.Tests.FEM;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private const int subdomainID = 0;

        //static void Main(string[] args)
        //{
        //    for (int example = 28; example < 29; example++)
        //    {
        //        CnstValues.exampleNo = example;
        //        CnstValues.runOnlyHexaModel = false;

        //        (Model model1, double[] uc1, Vector globalU1, bool IsFetiDpSolver3d) = SeparateCodeCheckingClass5b_bNEW_debugGit.RunExample();
        //        (Model model2, double[] uc2, Vector globalU2) = SeparateCodeCheckingClass5b_bNEW_debugGit.RunExampleSerial();
        //        Vector globalU1_2 = ReorderDirectSolverSolutionIn_globalU1_format(globalU1, globalU2, model1, model2);
        //        var check = ((globalU1 - globalU1_2).Norm2()) / (globalU1.Norm2());
        //        var check2 = (globalU1 - globalU1_2);
        //        printGlobalSolutionStats(check, IsFetiDpSolver3d);





        //        //CnstValues.runOnlyHexaModel = true;
        //        //CnstValues.preventOutputFileWrite(); 

        //        //FetiDP3dSolverSerialTestsInput.TestSolutionGlobalDisplacements();

        //        //CnstValues.RestoreDefaultBoolValues();
        //    }

        //}

        public static void MainMPI(string[] args)
        {
            //ProfileFetiDPCantileverBeam2D.Run();

            var suite = new MpiTestSuite();
            //suite.AddFact(materialParrallelExecutionTest.TestMaterialUpdateOnly);
            //suite.AddFact(materialManagerParrallelExecutionTest1.TestMaterialUpdateOnly);
            //suite.AddFact(Hexa8NonLinearCantileverDefGradDevelop4.ParallelNonLinearCantilever);
            //suite.AddFact(Hexa8NonLinearCantileverDefGradDevelop4Multiscale.ParallelNonLinearCantilever);
            suite.AddFact(Hexa8NonLinearCantileverDefGradDevelop4MultiscaleGraphene.ParallelNonLinearCantilever);
            suite.RunTests(args);
        }

        public static void Main(string[] args)
        {
            AssemblyCheck.WritingPhase = true;
            (var sttressesFeti1, var stressesFeti) =MemoryProfilerRveExample.CheckExample46InputInCodeRAMconsumption();
            AssemblyCheck.WritingPhase = false;
            (var sttressesFeti12, var stressesFeti2) =MemoryProfilerRveExample.CheckExample46InputInCodeRAMconsumptionV2();
            AssemblyCheck.ViewResults();
        }
        

        private static void printGlobalSolutionStats(double check, bool IsFetiDpSolver3d)
        {
            var cnstVal = new CnstValues();
            if (CnstValues.printGlobalSolutionStats)
            {
                string[] statsLines = new string[] { "GlobalSolutionError=" + check.ToString() + ",", };
                if (IsFetiDpSolver3d)
                {
                    var statsOutputPath = cnstVal.interfaceSolverStatsPath + @"\GlobalSolution_FetiDP_3d_stats.txt";
                    cnstVal.WriteToFileStringArray(statsLines, statsOutputPath);
                }
                else
                {
                    var statsOutputPath = cnstVal.interfaceSolverStatsPath + @"\GlobalSolution_FetiDP_stats.txt";
                    cnstVal.WriteToFileStringArray(statsLines, statsOutputPath);
                }
            }
        }

        private static Vector ReorderDirectSolverSolutionIn_globalU1_format(Vector globalU1, Vector globalU2 , Model model1, Model model2)
        {
            //var freeNodesIds = model1.NodesDictionary.Values.Where(x => x.Constraints.Count == 0).Select(x=>x.ID).ToList();
            var freeNodesIds = model1.GlobalDofOrdering.GlobalFreeDofs.GetRows().Select(x => x.ID);
            Vector globalU1_2 = Vector.CreateZero(globalU1.Length);
            foreach (int nodeID in freeNodesIds)
            {
                foreach (var dof in model2.GlobalDofOrdering.GlobalFreeDofs.GetDataOfRow(model2.GetNode(nodeID)).Keys)
                {
                    int model1Dof_order = model1.GlobalDofOrdering.GlobalFreeDofs[model1.GetNode(nodeID), dof];
                    int model2Dof_order = model2.GlobalDofOrdering.GlobalFreeDofs[model2.GetNode(nodeID), dof];

                    globalU1_2[model1Dof_order] = globalU2[model2Dof_order];
                }
            }

            return globalU1_2; 
        }

        private static void SolveBuildingInNoSoilSmall()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.NodesDictionary[21], DOF = StructuralDof.TranslationX });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.01, 0.1);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");

            //TODO: No loads have been defined so the result is bound to be 0.
        }

        private static void SolveCantileverWithStochasticMaterial()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }
    }
}
