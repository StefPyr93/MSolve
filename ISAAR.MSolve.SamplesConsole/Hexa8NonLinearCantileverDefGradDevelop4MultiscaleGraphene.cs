
using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Analyzers.ObjectManagers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MSAnalysis.remoteMatImplementations;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.SamplesConsole;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class Hexa8NonLinearCantileverDefGradDevelop4MultiscaleGraphene
    {
        private const int subdomainID = 1;
                    
        public static void ParallelNonLinearCantilever(int numProcesses)
        {
                   
            int numSubdomains = numProcesses;
            var procs = ProcessDistribution.CreateDistribution(numProcesses, numSubdomains); // FetiDPDofSeparatorMpiTests .CreateModelAndDofSeparator

            Console.Write("thread sleeping for sychronization");
            //Console.Write($"waiting time = " + procs.OwnRank*20000);
            System.Threading.Thread.Sleep(/*(procs.OwnRank+1)**/20000);

            #region rve material choice and manager creation and if the micro model will be solved serially or feti
            //CnstValues.exampleNo = 28;
            CnstValues.exampleNo = 46; CnstValues.parameterSet = ParameterSet.stiffLargerRve;
            CnstValues.runOnlyHexaModel = false;
            CnstValues.isInputInCode_forRVE = true;
            CnstValues.useInput_forRVE = true;
            CnstValues.PreventMATLABandTotalOutput();

            #region discretization data
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

            var rveBuilderSuitesparse = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number, false);
            var material1 = new MicrostructureDefGrad3D(rveBuilderSuitesparse,
                skylinemodel => (new SuiteSparseSolver.Builder()).BuildSolver(skylinemodel), false, 1);
            IMaterialManager materialManager = new MaterialManagerMpi2(material1, procs);

            //var rveBuilderFeti = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, true, mpgp,
            //subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);
            //var material1 = new MicrostructureDefGrad3DSerial(rveBuilderFeti,
            //    rveBuilderFeti.GetAppropriateSolverMpi, false, 1, false, false);
            //IMaterialManager materialManager = new MaterialManagerMpi2(material1, procs);

            #endregion

            Model model = new Model();//'

            if(procs.IsMasterProcess)
            {
                BuildCantileverModel(model, materialManager);
            }


            var increments = 2;//.
            var resTol = 1E-3;// sto multiscale xrhsimopoihsame thn default.
            int maxIters = 100;
            int itersRebuild = 1;


            (StaticAnalyzerDevelopMpi parentAnalyzer, LoadControlAnalyzerDevelop4Mpi childAnalyzer) = GetAnalyzers(procs, model, increments, maxIters,
                itersRebuild, resTol, materialManager);
            


            materialManager.Initialize();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            if (procs.IsMasterProcess)
            {
                
                IncrementalDisplacementsLog log1 = childAnalyzer.IncrementalDisplacementsLog;//
                IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
                bool isProblemSolvedCorrectly = AreDisplacementsSame(expectedDisplacements, log1,1E-6);
                if (isProblemSolvedCorrectly)
                {
                    Console.WriteLine($"Problem is solved correctly ");
                }
                else
                {
                    Console.WriteLine($"the problem has not been solved correctly");
                }
                //if (CnstValues.writeFe2MacroscaleSolution)
                //{
                //    double[][] solutionVectors = ExtractCalculatedSolutionsMacro(log1, subdomainID);
                //    for (int i1 = 0; i1 < solutionVectors.Length; i1++)
                //    {
                //        DdmCalculationsGeneral.WriteToFileVector(solutionVectors[i1], (new CnstValues()).exampleOutputPathGen + (@"\Msolve_solution\MacroscaleSolutionatInctrment" + i1 + ".txt"));
                //    }
                //}

            }


        }

        private static double[][] ExtractCalculatedSolutionsMacro(IncrementalDisplacementsLog log1, int subdomainId)
        {
            double[][] incrementSolutions = new double[log1.dofDisplacementsPerIter.Count][];
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < log1.dofDisplacementsPerIter.Count; iter++)
            {
                incrementSolutions[iter] = new double[log1.dofDisplacementsPerIter[iter][subdomainId].Keys.Count];
                int thesi = 0;
                foreach (int dof in log1.dofDisplacementsPerIter[iter][subdomainId].Keys)
                {
                    incrementSolutions[iter][thesi] = log1.dofDisplacementsPerIter[iter][subdomainId][dof];
                    thesi++;
                }
            }
            return incrementSolutions;
        }

        private static (StaticAnalyzerDevelopMpi parentAnalyzer, LoadControlAnalyzerDevelop4Mpi childAnalyzer) GetAnalyzers(ProcessDistribution procs,
            Model model, int increments, int maxIters, int itersRebuild, double resTol, IMaterialManager materialManager)
        {
            if (procs.IsMasterProcess)
            {
                // Solver
                var solverBuilder = new SkylineSolver.Builder();
                ISolver solver = solverBuilder.BuildSolver(model);
                // Problem type
                var provider = new ProblemStructural(model, solver);
                var childAnalyzer = new LoadControlAnalyzerDevelop4Mpi(model, solver, provider, increments, maxIters, itersRebuild, resTol, materialManager, procs);
                var parentAnalyzer = new StaticAnalyzerDevelopMpi(model, solver, provider, childAnalyzer, procs);
                var watchDofs = new Dictionary<int, int[]>();
                watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
                //var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
                //childAnalyzer.TotalDisplacementsPerIterationLog = log1; //.
                var log1 = new IncrementalDisplacementsLog(watchDofs);
                childAnalyzer.IncrementalDisplacementsLog = log1;

                return (parentAnalyzer, childAnalyzer);
            }
            else
            {
                var childAnalyzer = new LoadControlAnalyzerDevelop4Mpi(materialManager, procs, increments, maxIters, itersRebuild);
                var parentAnalyzer = new StaticAnalyzerDevelopMpi(childAnalyzer, procs);

                return (parentAnalyzer, childAnalyzer);
            }
        }

        public static void BuildCantileverModel(Model model, IMaterialManager materialManager)
        {
            double load_value = 0.00219881744271988174427;
            int subdomainID = 1;
            var subd1 = new Subdomain(subdomainID);
            model.SubdomainsDictionary.Add(subdomainID, subd1);
            IContinuumMaterial3DDefGrad material1 = new RemoteMaterial(materialManager);

            

            double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
            {0.250000,-0.250000,-1.000000},
            {-0.250000,0.250000,-1.000000},
            {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000},
            {0.250000,-0.250000,-0.500000},
            {-0.250000,0.250000,-0.500000},
            {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000},
            {0.250000,-0.250000,0.000000},
            {-0.250000,0.250000,0.000000},
            {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000},
            {0.250000,-0.250000,0.500000},
            {-0.250000,0.250000,0.500000},
            {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000},
            {0.250000,-0.250000,1.000000},
            {-0.250000,0.250000,1.000000},
            {0.250000,0.250000,1.000000}};

            int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14}, };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }

            // orismos elements 
            Element e1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1.ID, e1);
            }

            // constraint vashh opou z=-1
            for (int k = 1; k < 5; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
            }

            // fortish korufhs
            Load load1;
            for (int k = 17; k < 21; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = StructuralDof.TranslationX,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
        

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, TotalDisplacementsPerIterationLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, IncrementalDisplacementsLog computedDisplacements,double tol =1E-13)
        {
            var comparer = new ValueComparer(tol);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }


        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[2]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

            expectedDisplacements[0] = new Dictionary<int, double> {  {0, 0.032945773642662511},
                    {11, -0.026568607197976647 }, {23, -0.053540818543937239 }, {35, -0.077036613067688109}, {47, -0.09678117550794936}};


            expectedDisplacements[1] = new Dictionary<int, double> { {0, 0.065619830559031339},
                {11, -0.053046604784112175}, {23, -0.11346801809331361}, {35, -0.17571559014925475}, {47, -0.23723117394203658},};


            return expectedDisplacements;
        }



    }

}
