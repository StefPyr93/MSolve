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
//using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.SamplesConsole.SupportiveClasses;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.MSAnalysis.RveTemplates.SupportiveClasses;
using Xunit;

namespace ISAAR.MSolve.Tests.FEMpartB.SeparationBenchmarks2
{
    public static class SeparateCodeCheckingClass5b_bNEW_debugGitTest
    {//check unpushed changes commit

        [Fact]
        static void Solve()
        {
            for (int example = 28; example < 29; example++)
            {
                CnstValues.exampleNo = example;
                CnstValues.runOnlyHexaModel = false;

                (Model model1, double[] uc1, Vector globalU1, bool IsFetiDpSolver3d) = SeparateCodeCheckingClass5b_bNEW_debugGitTest.RunExample();
                (Model model2, double[] uc2, Vector globalU2) = SeparateCodeCheckingClass5b_bNEW_debugGitTest.RunExampleSerial();
                Vector globalU1_2 = ReorderDirectSolverSolutionIn_globalU1_format(globalU1, globalU2, model1, model2);
                var check = ((globalU1 - globalU1_2).Norm2()) / (globalU1.Norm2());
                var check2 = (globalU1 - globalU1_2);
                printGlobalSolutionStats(check, IsFetiDpSolver3d);


                double[] expectedGlobalSolutionFeti3Dfirst11values = new double[11] {0.000126146992609552, 0.000130114215483872, 8.33384821217284E-05,9.32792268779396E-05,9.74136328421431E-05,
                    5.26728981221024E-05,3.18832581261157E-05,3.26803860323064E-05,6.36506185081758E-05,-1.95427905842077E-05,-2.34221379496245E-05 };
                double[] calculatedGlobalSolutionFeti = (new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }).Select(x => globalU1[x]).ToArray();
                    


                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(expectedGlobalSolutionFeti3Dfirst11values, calculatedGlobalSolutionFeti));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(new double[] { check }, new double[] { 5.7123480244662094E-05 },0.2)); // actual error = 5.7123480244662094E-05



                #region old comments for hexa only model
                //CnstValues.runOnlyHexaModel = true;
                //CnstValues.preventOutputFileWrite(); 

                //FetiDP3dSolverSerialTestsInput.TestSolutionGlobalDisplacements();

                //CnstValues.RestoreDefaultBoolValues();
                #endregion
            }

        }

        [Fact]
        static void SolveInput()
        {
            for (int example = 28; example < 29; example++)
            {
                CnstValues.exampleNo = example;
                CnstValues.runOnlyHexaModel = false;
                CnstValues.isInputInCode_forRVE = true;
                CnstValues.PreventMATLABandTotalOutput();

                (Model model1, double[] uc1, Vector globalU1, bool IsFetiDpSolver3d) = SeparateCodeCheckingClass5b_bNEW_debugGitTest.RunExample();
                (Model model2, double[] uc2, Vector globalU2) = SeparateCodeCheckingClass5b_bNEW_debugGitTest.RunExampleSerial();
                Vector globalU1_2 = ReorderDirectSolverSolutionIn_globalU1_format(globalU1, globalU2, model1, model2);
                var check = ((globalU1 - globalU1_2).Norm2()) / (globalU1.Norm2());
                var check2 = (globalU1 - globalU1_2);
                printGlobalSolutionStats(check, IsFetiDpSolver3d);


                double[] expectedGlobalSolutionFeti3Dfirst11values = new double[11] {0.000126146992609552, 0.000130114215483872, 8.33384821217284E-05,9.32792268779396E-05,9.74136328421431E-05,
                    5.26728981221024E-05,3.18832581261157E-05,3.26803860323064E-05,6.36506185081758E-05,-1.95427905842077E-05,-2.34221379496245E-05 };
                double[] calculatedGlobalSolutionFeti = (new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }).Select(x => globalU1[x]).ToArray();



                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(expectedGlobalSolutionFeti3Dfirst11values, calculatedGlobalSolutionFeti,1e-8));
                Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(new double[] { check }, new double[] { 5.7123480488295326E-05 }, 0.2)); // error 5.7123480488295326E-05



                #region old comments for hexa only model
                //CnstValues.runOnlyHexaModel = true;
                //CnstValues.preventOutputFileWrite(); 

                //FetiDP3dSolverSerialTestsInput.TestSolutionGlobalDisplacements();

                //CnstValues.RestoreDefaultBoolValues();
                #endregion
            }

        }



        //prosthiki model.ConnectDataStructures entos rve gia na vrei to output node.Subdomains =/=0
        public static (Model, double[], Vector, bool) RunExample()
        {
            // EPILOGH RVE
            int subdiscr1;//= 4;// 4;// 6;
            int discr1;//= 2;// 3;//4;

            int discr3;//= discr1 * subdiscr1;// 23;
            int subdiscr1_shell;//= 6;//14;
            int discr1_shell;// = 1;
            int graphene_sheets_number;// =2; //periektikothta 0.525% 
            double scale_factor;//= 1; //PROSOXH

            (subdiscr1, discr1, subdiscr1_shell, discr1_shell, graphene_sheets_number, scale_factor) = GetGrRveExampleDiscrDataFromFile(new CnstValues());
            discr3 = discr1 * subdiscr1;

            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;

            bool run_new_corner = true;
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3D(1, true, mpgp,
            //subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);
            CnstValues.useInput_forRVE = true;
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            //rveBuilder.useInput = true;
            // EPILOGH RVE
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateDevelop(1, true); //edw ginetai develop h feti dp gia provlhmata 3d
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateLARGE(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop(1, true);

            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelopHSTAM(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopb(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbLARGE(1, true);
            //var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elem(1, true);

            bool WRITESTIFFNESSES = CnstValues.WRITESTIFFNESSES;

            #region model nodes and load
            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            //model.ConnectDataStructures();

            double load_value = 1; //A.7
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);
            var cornerNodesAndSubds = rveBuilder.CornerNodesIdAndsubdomains;
            model.ConnectDataStructures();

            Dictionary<int, HashSet<INode>> cornerNodes = rveBuilder.cornerNodes;
            Dictionary<ISubdomain, HashSet<INode>> extraConstrNodesofsubd = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subd in model.EnumerateSubdomains()) extraConstrNodesofsubd.Add((ISubdomain)subd, new HashSet<INode>());
            foreach(var extraConstrNodeList in rveBuilder.extraConstraintsNoeds)
            {
                int extraNodeId = extraConstrNodeList[0];
                foreach (var subd in model.NodesDictionary[extraNodeId].SubdomainsDictionary.Values)
                {
                    extraConstrNodesofsubd[subd].Add(model.NodesDictionary[extraNodeId]);
                }
            }
            if(WRITESTIFFNESSES)
            {
                var subdExtraConstrsFirstNode = extraConstrNodesofsubd.Select(x => new KeyValuePair<int, int[]>  ( x.Key.ID, x.Value.Select(y => y.ID).ToArray())  ).ToDictionary(x=>x.Key,x=>x.Value); 
                DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdExtraConstrsFirstNode, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\subdExtraConstrsFirstNode.txt");
            }
            #endregion



            #region setup solver problem and initialize
            //Setup solver
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = 1E-4,
                MaxIterationsProvider = new FixedMaxIterationsProvider(1000)
            };
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();

            var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x => x.Item1, x => x.Value);

            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);
            var midSideNodeSelection = new UserDefinedMidsideNodes(extraConstrNodesofsubd, new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
            
            //var fetiSolverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);  //A.3
            //var matrixManagerFactory = new FetiDP3dMatrixManagerFactoryDense();   //A.3.1
            //var matrixManagerFactory = new FetiDP3dMatrixManagerFactorySkyline();   //A.3.1
            var matrixManagerFactory = new FetiDP3dMatrixManagerFactorySuiteSparse(new OrderingAmdSuiteSparse());   //A.3.1
            var fetiSolverBuilder = new FetiDP3dSolverSerial.Builder(matrixManagerFactory);  //A.3

            //fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.StiffnessDistribution = StiffnessDistributionType.HeterogeneousCondensed;
            fetiSolverBuilder.Preconditioning = new DirichletPreconditioning();
            fetiSolverBuilder.PcgSettings = pcgSettings;

            Crosspoints crosspoints = Crosspoints.FullyRedundant; //A.2
            //Crosspoints crosspoints = Crosspoints.Minimum; //A.2

            ICrosspointStrategy crosspointStrategy;
            crosspointStrategy = new MinimumConstraints();
            if (crosspoints == Crosspoints.FullyRedundant) crosspointStrategy = new FullyRedundantConstraints();
            fetiSolverBuilder.CrosspointStrategy = crosspointStrategy;

            //FetiDPSolverSerial fetiSolver = fetiSolverBuilder.Build(model, cornerNodeSelection); //A.1
            FetiDP3dSolverSerial fetiSolver = fetiSolverBuilder.Build(model, cornerNodeSelection, midSideNodeSelection); //A.1

            //FetiDPSolverPrint fetiSolver = fetiSolverBuilder.BuildSolver(model);
            //model.ConnectDataStructures();

            // Run the analysis
            //var problem = new ProblemStructural(model, (ISolver)fetiSolver);
            //var linearAnalyzer = new LinearAnalyzer(model, (ISolver)fetiSolver, problem);
            //var staticAnalyzer = new StaticAnalyzer(model, (ISolver)fetiSolver, problem, linearAnalyzer);
            //staticAnalyzer.Initialize();
            RunAnalysis(model, fetiSolver); // apo to paradeigma ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests.FetiDPSolverSerialTests.RunAnalysis()
            #endregion


            #region print subdomain stiffness matrix and global dofs
            // @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}Iter{1}Stiffness.txt";
            string print_path_gen = rveBuilder.subdomainOutputPath + @"\subdomain_matrices_and_data\Subdomain{0}Iter{1}Stiffness.txt";
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                //var subdMatrix= provider.CalculateMatrix(subdomain);

                string subdID = subdomain.ID.ToString();
                var subdMatrix = fetiSolver.GetLinearSystem(subdomain).Matrix; // subdomainMatrixes[subdomain.ID];

                string counter_data = 1.ToString();
                string print_path = string.Format(print_path_gen, subdID, counter_data);

                var writer = new MatlabWriter();
                if (WRITESTIFFNESSES) writer.WriteToFile((ISparseMatrix)subdMatrix, print_path, false);
            }
            //string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}GlobalDofs.txt";
            string print_path_gen2 = rveBuilder.subdomainOutputPath + @"\subdomain_matrices_and_data\Subdomain{0}GlobalDofs.txt";
            foreach (Subdomain subdomain in model.EnumerateSubdomains())
            {
                double[] subdomainGlobalDofs = new double[subdomain.FreeDofOrdering.NumFreeDofs];

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };

                foreach (Node node in subdomain.Nodes.Values)
                {
                    foreach (StructuralDof dof in dofs)
                    {
                        bool check = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(node, dof, out int dofValue);
                        if (check)
                        {
                            subdomainGlobalDofs[dofValue] = model.GlobalDofOrdering.GlobalFreeDofs[node, dof];

                            if (model.GlobalDofOrdering.GlobalFreeDofs[node, dof] == 0)
                            {
                                string breakpoint = "here";
                            }
                        }
                    }
                }

                string subdID = subdomain.ID.ToString();
                string print_path = string.Format(print_path_gen2, subdID);
                var writer = new MatlabWriter();
                if (WRITESTIFFNESSES) writer.WriteToFile(Vector.CreateFromArray(subdomainGlobalDofs), print_path, false);
            }
            #endregion

            #region print: corner nodes and dofIds BR nodes and dofIds
            Dictionary<int, int[]> CornerNodesIdAndGlobalDofs = new Dictionary<int, int[]>(rveBuilder.CornerNodesIds.Keys.Count());//nodeID, globalDofs
            Dictionary<int, int[]> subdBRNodesAndGlobalDOfs = new Dictionary<int, int[]>(rveBuilder.subdFreeBRNodes.Keys.Count());//nodeID, globalDofs
            foreach (int corrnerNodeID in rveBuilder.CornerNodesIds.Keys)
            {
                Node CornerNode = model.NodesDictionary[corrnerNodeID];
                CornerNodesIdAndGlobalDofs.Add(corrnerNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationZ]});

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(CornerNode, StructuralDof.RotationX, out int globalDofId3);
                if (check)
                {
                    string breakpoint = "here";
                }
            }
            foreach (int boundaryNodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node boundaryNode = model.NodesDictionary[boundaryNodeID];

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(boundaryNode, StructuralDof.RotationX, out int globalDofId4);
                if (!check)
                {
                    subdBRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ]});
                }
                else
                {
                    subdBRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[5] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationY]});
                }
            }
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndGlobalDofs, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\CornerNodesAndGlobalDofsIds.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdBRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\SubdBRNodesAndGlobalDofsIds.txt");
            #endregion

            #region overwrite extra constraint nodes order and print them
            List<List<int>> extraConstraintsNoeds = rveBuilder.extraConstraintsNoeds;
            bool changeOrder = false;
            if (changeOrder)
            {
                extraConstraintsNoeds[0] = new List<int>() { 38 };
                extraConstraintsNoeds[1] = new List<int>() { 58 };
                extraConstraintsNoeds[2] = new List<int>() { 62 };
                extraConstraintsNoeds[3] = new List<int>() { 64 };
                extraConstraintsNoeds[4] = new List<int>() { 70 };
                extraConstraintsNoeds[5] = new List<int>() { 100 };
            }
            

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = GetExtraConstrNodesPositions(subdBRNodesAndGlobalDOfs, extraConstraintsNoeds, model);
            Dictionary<int, int[]> ExtraConstrAveIdsAndTheirBRNodesTheseis = GetExtraConstrNodesPositions(subdBRNodesAndGlobalDOfs, rveBuilder.extraConstraintsNoedsAve, model);

            bool ommitZeros = true;
            if (ommitZeros) { ExtraConstrIdAndTheirBRNodesTheseis = OmmitZeros(ExtraConstrIdAndTheirBRNodesTheseis); }
            if (ommitZeros) { ExtraConstrAveIdsAndTheirBRNodesTheseis = OmmitZeros(ExtraConstrAveIdsAndTheirBRNodesTheseis); }

            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBRNodesTheseis, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\ExtraConstrIdAndTheirBRNodesTheseis.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrAveIdsAndTheirBRNodesTheseis, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\ExtraConstrAveIdsAndTheirBRNodesTheseis.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.printNodeData(model, extraConstraintsNoeds, rveBuilder.CornerNodesIdAndsubdomains, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\ExtraNodeCoordinates.txt", @"\subdomain_matrices_and_data\CornerNodeCoordinates.txt");
            #endregion
                                   
            #region coupled data arrays
            Dictionary<int, int[]> GlobalDofCoupledDataSubdIds = new Dictionary<int, int[]>(3 * (CornerNodesIdAndGlobalDofs.Count() + subdBRNodesAndGlobalDOfs.Count));
            Dictionary<int, int[]> GlobalDofCoupledDataLocalDofsInSubdIds = new Dictionary<int, int[]>(3 * (CornerNodesIdAndGlobalDofs.Count() + subdBRNodesAndGlobalDOfs.Count));
            foreach (int nodeID in rveBuilder.CornerNodesIds.Keys)
            {
                Node node = model.NodesDictionary[nodeID];
                int[] subdIds = node.SubdomainsDictionary.Keys.ToArray();

                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                foreach (StructuralDof doftype in dofs)
                {
                    int globalDofId1 = model.GlobalDofOrdering.GlobalFreeDofs[node, doftype];
                    int[] localIds = new int[subdIds.Length];

                    for (int i1 = 0; i1 < subdIds.Length; i1++)
                    {
                        localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                    }
                    GlobalDofCoupledDataSubdIds.Add(globalDofId1, subdIds);
                    GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId1, localIds);
                }
            }
            foreach (int nodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node node = model.NodesDictionary[nodeID];
                int[] subdIds = node.SubdomainsDictionary.Keys.ToArray();

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };

                foreach (StructuralDof doftype in dofs)
                {
                    bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, doftype, out int globalDofId2);
                    if (check)
                    {
                        int[] localIds = new int[subdIds.Length];

                        for (int i1 = 0; i1 < subdIds.Length; i1++)
                        {
                            localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                        }
                        GlobalDofCoupledDataSubdIds.Add(globalDofId2, subdIds);
                        GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId2, localIds);
                    }
                }

            }

            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataSubdIds, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\GlobalDofCoupledDataSubdIds.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataLocalDofsInSubdIds, rveBuilder.subdomainOutputPath, @"\subdomain_matrices_and_data\GlobalDofCoupledDataLocalDofsInSubdIds.txt");
            #endregion
            

            //staticAnalyzer.Solve();


            #region  Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in model.EnumerateSubdomains()) sudomainDisplacements[ls.ID] = fetiSolver.GetLinearSystem(ls).Solution;
            Vector globalU = fetiSolver.GatherGlobalDisplacements();// sudomainDisplacements);

            Node monitoredNode = model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key];
            int globalDofId = model.GlobalDofOrdering.GlobalFreeDofs[monitoredNode, StructuralDof.TranslationZ];
            double solution = globalU[globalDofId];


            double[] uc = new double[3 * cornerNodesAndSubds.Count()];

            int node_counter = 0;
            foreach (int nodeId in cornerNodesAndSubds.Keys)
            {
                //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                Node node = model.NodesDictionary[nodeId];
                for (int i1 = 0; i1 < 3; i1++)
                {
                    int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
                    uc[3 * node_counter + i1] = globalU[globalDof];

                }
                node_counter++;
            }
            //(new ISAAR.MSolve.LinearAlgebra.Output.Array1DWriter()).WriteToFile(globalU.CopyToArray(), rveBuilder.subdomainOutputPath + @"\Msolve_solution\Global_solution.txt");
            bool IsFetiDpSolver3d = false;
            bool writeFetiSolutionVectors = CnstValues.writeFetiSolutionVectors;
            if (writeFetiSolutionVectors)
            {
                if (fetiSolver is FetiDP3dSolverSerial)
                {
                    DdmCalculationsGeneral.WriteToFileVector(globalU.CopyToArray(), rveBuilder.subdomainOutputPath + @"\Msolve_solution\Global_solution_fetiDP3D.txt");
                    DdmCalculationsGeneral.WriteToFileVector(uc, rveBuilder.subdomainOutputPath + @"\Msolve_solution\Corner_solution_fetiDP3D.txt");
                    IsFetiDpSolver3d = true;
                }
                else
                {
                    DdmCalculationsGeneral.WriteToFileVector(globalU.CopyToArray(), rveBuilder.subdomainOutputPath + @"\Msolve_solution\Global_solution_fetiDP.txt");
                    DdmCalculationsGeneral.WriteToFileVector(uc, rveBuilder.subdomainOutputPath + @"\Msolve_solution\Corner_solution_fetiDP.txt");
                    IsFetiDpSolver3d = true;
                }
            }
            #endregion

            #region  overwrite data model region
            bool run_overwrite_data_region = CnstValues.run_overwrite_data_region;
            bool print_hexa_model = CnstValues.print_hexa_model;
            if (run_overwrite_data_region)
            {
                if (print_hexa_model)
                { PrintHexaModelData(model, rveBuilder.subdomainOutputPath); }

                Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesIds = GetExtraConstrNodesIds(subdBRNodesAndGlobalDOfs, extraConstraintsNoeds, model);
                int[] brNodesMsolveWise = PrintUtilities.ReadIntVector(rveBuilder.subdomainOutputPath + @"\model_overwrite\subdomain_data_solver\RB_Nodes_IDs_MSOLVE_wise" + ".txt");
                Dictionary<int, int[]> subdBRNodesMsolveWiseAndGlobalDOfs = new Dictionary<int, int[]>();
                for (int i1 = 0; i1 < brNodesMsolveWise.Length; i1++) subdBRNodesMsolveWiseAndGlobalDOfs.Add(brNodesMsolveWise[i1], new int[1]);
                Dictionary<int, int[]> ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis = GetExtraConstrNodesPositions(subdBRNodesMsolveWiseAndGlobalDOfs, extraConstraintsNoeds, model);


                Dictionary<int, int[]> subd_msolve_BRNodesAndGlobalDOfs = new Dictionary<int, int[]>(rveBuilder.subdFreeBRNodes.Keys.Count());
                foreach (int boundaryNodeID in subdBRNodesMsolveWiseAndGlobalDOfs.Keys)
                {
                    Node boundaryNode = model.NodesDictionary[boundaryNodeID];

                    bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(boundaryNode, StructuralDof.RotationX, out int globalDofId4);
                    if (!check)
                    {
                        subd_msolve_BRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ]});
                    }
                    else
                    {
                        subd_msolve_BRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[5] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationY]});
                    }
                }
                if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subd_msolve_BRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\model_overwrite\subd_msolve_BRNodesAndGlobalDOfs.txt");
                if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBRNodesIds, rveBuilder.subdomainOutputPath, @"\model_overwrite\ExtraConstrIdAndTheirBRNodesIds.txt");
                if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis, rveBuilder.subdomainOutputPath, @"\model_overwrite\ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis.txt");
            }
            #endregion

            return (model, uc, globalU, IsFetiDpSolver3d);

        }

        private static (int subdiscr1, int discr1, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number, double scale_factor)  GetGrRveExampleDiscrDataFromFile(CnstValues cnstValues)
        {
            //DUPLICATE CHANGES IN SAMPLE CONSOLE
            if (!CnstValues.isInputInCode_forRVE)
            {
                int[] discrData = PrintUtilities.ReadIntVector(cnstValues.exampleDiscrInputPathGen + @"\subdiscr1_discr1_ subdiscr1_shell_discr1_shell_graphene_sheets_number" + ".txt");
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

        internal static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis
            solver.OrderDofs(false);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
            solver.BuildGlobalMatrix(new ElementStructuralStiffnessProvider());
            model.ApplyLoads();
            LoadingUtilities.ApplyNodalLoads(model, solver);
            solver.Solve();
        }

        private static Dictionary<int, int[]> OmmitZeros(Dictionary<int, int[]> extraConstrIdAndTheirBRNodesTheseis)
        {
            Dictionary<int, int[]> updaatedList = new Dictionary<int, int[]>();
            int counter = 0;
            foreach(KeyValuePair<int,int[]> extraNodesData in extraConstrIdAndTheirBRNodesTheseis)
            {
                int id = extraNodesData.Key;
                var initialTheseis = extraNodesData.Value.ToList();
                initialTheseis.RemoveAll(x => x == 0);
                var updated = initialTheseis.ToArray();
                updaatedList[id] = updated;
            }
            return updaatedList;
        }

        private static void PrintHexaModelData(Model model, string subdomainOutputPath)
        {
            Dictionary<int, List<int>> ElementIdsAndModelIds = new Dictionary<int, List<int>>();
            foreach (var element in model.EnumerateElements().Where(x=>(x.ElementType is Hexa8NonLinear| x.ElementType is Hexa8NonLinearDefGrad)))
            {
                List<int> ElementNodesIds = element.Nodes.Select(x => x.ID).ToList();
                ElementIdsAndModelIds.Add(element.ID, ElementNodesIds);
            }

            int[] ElementIds = model.EnumerateElements().Where(x => (x.ElementType is Hexa8NonLinear | x.ElementType is Hexa8NonLinearDefGrad)).Select(x => x.ID).ToArray();
            int[] subdomainIds = model.EnumerateSubdomains().Select(x => x.ID).ToArray();
            

            List<int> nodeIds = new List<int>();
            foreach (var hexaElementId in ElementIds)
            {
                nodeIds = nodeIds.Union(model.GetElement(hexaElementId).Nodes.Select(x => x.ID).ToList()).ToList();
            }
            int[] NodeIds = nodeIds.ToArray();

            int[,] ElementNodes = new int[ElementIds.GetLength(0), 8];
            int thesi = 0;
            foreach (var elementID in ElementIds)
            {
                int thesi2 = 0;
                foreach(Node node in model.ElementsDictionary[elementID].Nodes)
                {
                    ElementNodes[thesi, thesi2] = node.ID;
                    thesi2++;
                }
                thesi++;
            }

            double[,] NodeCoordinates = new double[model.NodesDictionary.Values.Count, 3];
            thesi = 0;
            foreach (var nodeID in NodeIds)
            {
                NodeCoordinates[thesi, 0] = model.NodesDictionary[nodeID].X1;
                NodeCoordinates[thesi, 1] = model.NodesDictionary[nodeID].X2;
                NodeCoordinates[thesi, 2] = model.NodesDictionary[nodeID].X3;
                thesi++;
            }

            //int[,] SubdElements = new int[model.EnumerateSubdomains().Count(), model.EnumerateSubdomains().ElementAt(0).EnumerateElements().Count()];
            Dictionary<int, int[]> subdElements = new Dictionary<int, int[]>(model.EnumerateSubdomains().Count());
            thesi = 0;
            foreach (var SubdId in subdomainIds)
            {
                int[] elements = model.SubdomainsDictionary[SubdId].Elements.Values.Where(x => (x.ElementType is Hexa8NonLinear | x.ElementType is Hexa8NonLinearDefGrad)).Select(x => x.ID).ToArray();
                subdElements.Add(SubdId, elements);
                //int thesi2 = 0;
                //foreach (var element in model.SubdomainsDictionary[SubdId].Elements.Values.Where(x => (x.ElementType is Hexa8NonLinear | x.ElementType is Hexa8NonLinearDefGrad)))
                //{
                //    var elementID = element.ID;
                //    SubdElements[thesi,thesi2]=elementID;
                //    thesi2++;
                //}
                //thesi++;

            }

            double[] elementStiffnessFactor = DetermineElementStiffnessFactor(model, ElementIds);

            List<int> constrainedNodes = new List<int>();
            foreach (var constraint in model.Constraints)
            {
                if (!(constrainedNodes.Contains(constraint.row.ID)))
                { constrainedNodes.Add(constraint.row.ID); }
            }
            var constraintIds = constrainedNodes.ToArray();


            bool outputTypeOriginal = false;

            if (outputTypeOriginal)
            {
                //print model reconstruction data 
                PrintUtilities.WriteToFileVectorMsolveInput(ElementIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementIds.txt");
                PrintUtilities.WriteToFileVectorMsolveInput(subdomainIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\subdomainIds.txt");
                PrintUtilities.WriteToFileVectorMsolveInput(NodeIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeIds.txt");
                PrintUtilities.WriteToFileVectorMsolveInput(constraintIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\constraintIds.txt");



                PrintUtilities.WriteToFileMsolveInput(NodeCoordinates, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeCoordinates.txt");
                PrintUtilities.WriteToFileMsolveInput(ElementNodes, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementNodes.txt");
                //PrintUtilities.WriteToFileMsolveInput(SubdElements, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\SubdElements.txt");
                PrintUtilities.WriteToFileDictionaryMsolveInput(subdElements, subdomainOutputPath, @"\model_overwrite\MsolveModel\subdElements.txt");
                PrintUtilities.WriteToFileVectorMsolveInput(elementStiffnessFactor, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementStiffnessFactors.txt");
            }
            else
            {
                //DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis, rveBuilder.subdomainOutputPath, @"\model_overwrite\ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis.txt");
                //int[] brNodesMsolveWise = ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(rveBuilder.subdomainOutputPath + @"\model_overwrite\subdomain_data_solver\RB_Nodes_IDs_MSOLVE_wise" + ".txt");
                //print model reconstruction data 
                PrintUtilities.WriteToFileVector(ElementIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementIds.txt");
                PrintUtilities.WriteToFileVector(subdomainIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\subdomainIds.txt");
                PrintUtilities.WriteToFileVector(NodeIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeIds.txt");
                PrintUtilities.WriteToFileVector(constraintIds, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\constraintIds.txt");



                PrintUtilities.WriteToFile(ElementNodes, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementNodes.txt");
                PrintUtilities.WriteToFile(NodeCoordinates, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeCoordinates.txt");
                //PrintUtilities.WriteToFileMsolveInput(SubdElements, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\SubdElements.txt");
                DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdElements, subdomainOutputPath, @"\model_overwrite\MsolveModel\subdElements.txt");
                PrintUtilities.WriteToFileVector(elementStiffnessFactor, subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementStiffnessFactors.txt");
            }



        }

        private static double[] DetermineElementStiffnessFactor(Model model, int[] elementIds)
        {
            var hexaElements = elementIds.Select(x => model.ElementsDictionary[x]).ToList();
            var nodes = hexaElements.Select(x => x.Nodes.OrderBy(y => AllignedDistance(y)).ElementAt(0)).ToArray();
            var uniqueNodes = nodes.Distinct();
            var uniqueNodesMultiplicities = new Dictionary<Node, int>(uniqueNodes.Count());

            foreach (var element in hexaElements)
            {
                var node = element.Nodes.OrderBy(y => AllignedDistance(y)).ElementAt(0);
                if (uniqueNodesMultiplicities.Keys.Contains(node))
                {
                    uniqueNodesMultiplicities[node]++;
                }
                else
                {
                    uniqueNodesMultiplicities.Add(node, 1);
                }


            }

            double[] elementStiffnessFactor = new double[elementIds.Length];
            for (int i1 = 0; i1 < hexaElements.Count(); i1++)
            {
                var node = nodes[i1];
                var multiplicitiy = uniqueNodesMultiplicities[node];
                elementStiffnessFactor[i1] = ((double)1 / multiplicitiy);
            }

            return elementStiffnessFactor;
        }

        private static double AllignedDistance(Node node)
        {
            //einai eswteriko ginomeno dld provolh sto dianusma (0,0,0)--->(1,1,1)
            return (node.X + node.Y + node.Z);  //Vector.CreateFromArray(new double[] { node.X, node.Y, node.Z }) * Vector.CreateFromArray(new double[] { 1, 1, 1 });
        }

        private static Dictionary<int, int[]> GetExtraConstrNodesPositions(Dictionary<int, int[]> subdBRNodesAndGlobalDOfs, List<List<int>> extraConstraintsNoedsIds, Model model)
        {
            int[] positionsOfBRNodes = new int[subdBRNodesAndGlobalDOfs.Keys.Max()];

            for (int i1 = 0; i1 < subdBRNodesAndGlobalDOfs.Keys.Count(); i1++)
            {
                int BRnodeID = subdBRNodesAndGlobalDOfs.ElementAt(i1).Key;
                positionsOfBRNodes[BRnodeID - 1] = i1 + 1;
            }

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = new Dictionary<int, int[]>();

            for (int i1 = 0; i1 < extraConstraintsNoedsIds.Count(); i1++)
            {
                int[] BRnodesTheseis = new int[extraConstraintsNoedsIds.ElementAt(i1).Count()];
                for (int i2 = 0; i2 < extraConstraintsNoedsIds.ElementAt(i1).Count(); i2++)
                {
                    BRnodesTheseis[i2] = positionsOfBRNodes[extraConstraintsNoedsIds.ElementAt(i1).ElementAt(i2) - 1];
                }

                ExtraConstrIdAndTheirBRNodesTheseis.Add(i1 + 1, BRnodesTheseis);
            }

            return ExtraConstrIdAndTheirBRNodesTheseis;
        }

        private static Dictionary<int, int[]> GetExtraConstrNodesIds(Dictionary<int, int[]> subdBRNodesAndGlobalDOfs, List<List<int>> extraConstraintsNoedsIds, Model model)
        {
            int[] positionsOfBRNodes = new int[subdBRNodesAndGlobalDOfs.Keys.Max()];

            for (int i1 = 0; i1 < subdBRNodesAndGlobalDOfs.Keys.Count(); i1++)
            {
                int BRnodeID = subdBRNodesAndGlobalDOfs.ElementAt(i1).Key;
                positionsOfBRNodes[BRnodeID - 1] = i1 + 1;
            }

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = new Dictionary<int, int[]>();

            for (int i1 = 0; i1 < extraConstraintsNoedsIds.Count(); i1++)
            {
                int[] BRnodesTheseis = new int[extraConstraintsNoedsIds.ElementAt(i1).Count()];
                for (int i2 = 0; i2 < extraConstraintsNoedsIds.ElementAt(i1).Count(); i2++)
                {
                    BRnodesTheseis[i2] = extraConstraintsNoedsIds.ElementAt(i1).ElementAt(i2);
                }

                ExtraConstrIdAndTheirBRNodesTheseis.Add(i1 + 1, BRnodesTheseis);
            }

            return ExtraConstrIdAndTheirBRNodesTheseis;
        }
        //needs to be corrected rve_multiple -> b kai to path kai ta stoixeia diakritopoihshs pou einai afhmena exwterika (Genika elegxoume connectDataStructures kai defineAppropriateConstraintsForBoundaryNodes)
        public static (Model, double[], Vector) RunExampleSerial()
        {
            // EPILOGH RVE
            int subdiscr1; //= 4;// 4;// 6;
            int discr1; //= 2;// 3;//4;
            // int discr2 dn xrhsimopoieitai
            int discr3;// = discr1 * subdiscr1;// 23;
            int subdiscr1_shell;// = 6;//14;
            int discr1_shell;// = 1;
            int graphene_sheets_number;// = 2; //periektikothta 0.525% 
            double scale_factor;// = 1;

            (subdiscr1, discr1, subdiscr1_shell, discr1_shell, graphene_sheets_number, scale_factor) = GetGrRveExampleDiscrDataFromFile(new CnstValues());
            discr3 = discr1 * subdiscr1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);

            (subdiscr1, discr1, subdiscr1_shell, discr1_shell, graphene_sheets_number, scale_factor) = GetGrRveExampleDiscrDataFromFile(new CnstValues());
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            if (CnstValues.parameterSet == ParameterSet.stiffCase)
            { mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90; }
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;

            CnstValues.useInput_forRVE = true;
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGitSerial(1, false, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number,false);
            //rveBuilder.useInput = true;

            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateLARGE(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopb(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbLARGE(1, false); // diorthose kai to parakatw path apla gia na mhn xtupaei.
            //var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elem(1, false); //A.1

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            var boundaryNodes = ModelAndNodes.Item2;

            var renumbering_vector_path = rveBuilder.renumbering_vector_path;

            
            //var mpgp = rveBuilder.mpgp; 
            var mp = mpgp.Item1; 
            var gp = mpgp.Item2;
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            Dictionary<int, double[]> CornerNodesIds;
            Dictionary<int, int[]> CornerNodesIdAndsubdomains;
            int[][] CornerNodesData = rveBuilder.CornerNodesData;
            CornerNodesIds = rveBuilder.CornerNodesIds;
            CornerNodesIdAndsubdomains = rveBuilder.CornerNodesIdAndsubdomains;

            //A.6
            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);

            //A.7
            bool boundaryNodesOrPaktwsh=true;
            if (boundaryNodesOrPaktwsh)
            {
                DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);
            }
            else
            {
                renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
                int[][] paktwsiNodesData = new int[4][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
                int thesi = 0;
                int j1 = 0;
                for (int i2 = 0; i2 < 2; i2++)
                {
                    for (int i3 = 0; i3 < 2; i3++)
                    {
                        paktwsiNodesData[thesi] = new int[3] { j1, i2 * 1, i3 * 1 };
                        thesi++;
                    }
                }
                Dictionary<int, Node> paktwshAristeraNodes = new Dictionary<int, Node>();
                for (int i1 = 0; i1 < paktwsiNodesData.Length; i1++)
                {
                    int h1 = paktwsiNodesData[i1][0]; int h2 = paktwsiNodesData[i1][1]; int h3 = paktwsiNodesData[i1][2];
                    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                    paktwshAristeraNodes.Add(nodeID, model.NodesDictionary[nodeID]);


                }
                DefineAppropriateConstraintsForBoundaryNodes(model, paktwshAristeraNodes);
            }


            #region define solver from Quad4LinearCantilever %81 and IntegrationElastic... %37 tests.
            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            // Request output
            //int[] outputPositions = new int[model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs];
            //for( int i1 = 0;  i1 < model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs ; i1++ )
            //{

            //}
            //childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { 0 });

            // Run the anlaysis - part a
            parentAnalyzer.Initialize();

            //A.8    
            bool print_data = false;
            if (print_data)
            {
                printElementStiffnessAndData(model);
            }

            // Run the anlaysis - part b
            parentAnalyzer.Solve();
            #endregion

            var globalU = solver.LinearSystems[1].Solution.CopyToArray();// fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
            
            double[] uc = new double[3 * CornerNodesIds.Count()];

            int node_counter = 0;
            foreach (int nodeId in CornerNodesIds.Keys)
            {
                //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                Node node = model.NodesDictionary[nodeId];
                for (int i1 = 0; i1 < 3; i1++)
                {
                    int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
                    uc[3 * node_counter + i1] = globalU[globalDof];

                }
                node_counter++;
            }
            if (CnstValues.writeFetiSolutionVectors)
            {
                DdmCalculationsGeneral.WriteToFileVector(globalU, rveBuilder.subdomainOutputPath + @"\Msolve_solution\Global_solution_Direct.txt");
                DdmCalculationsGeneral.WriteToFileVector(uc, rveBuilder.subdomainOutputPath + @"\Msolve_solution\Corner_solution_Direct.txt");
            }

            return (model, uc, Vector.CreateFromArray(globalU));
        }

        public static void printElementStiffnessAndData(Model model)
        {
            var dofOrdering = model.SubdomainsDictionary[0].FreeDofOrdering;
            var matrixProvider = new ElementStructuralStiffnessProvider();
            
            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}Stiffness.txt";
            string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}elementDofIndices.txt";
            string print_path_gen3 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}subdomainDofIndices.txt";
            foreach (Element element in model.ElementsDictionary.Values)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                var elementMatrix = matrixProvider.Matrix(element).CopytoArray2D();
                string counterData = element.ID.ToString();
                string print_path = string.Format(print_path_gen, counterData);
                var writer = new Array2DWriter();
                writer.WriteToFile(elementMatrix, print_path, false);

                //Aray1Dwriter
                string print_path2 = string.Format(print_path_gen2, counterData);
                var writer2 = new MatlabWriter();
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(elementDofIndices)), print_path2, false);

                string print_path3 = string.Format(print_path_gen3, counterData);
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(subdomainDofIndices)), print_path3, false);
            }
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            //// DUPLICATE ANNY CHANGES IN ONERVEEXAMPLEMPI
            //// kai sto integration Separatecodecheckingclass....GitTest

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

        private static void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        public static double[] Cnvrt(int[] array)
        {
            double[] array2 = new double[array.Length];
            for(int i1=0; i1<array.Length; i1++)
            {
                array2[i1] = (double)array[i1];
            }
            return array2;
        }

        public enum Crosspoints { Minimum, FullyRedundant }

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

        private static Vector ReorderDirectSolverSolutionIn_globalU1_format(Vector globalU1, Vector globalU2, Model model1, Model model2)
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


    }
}
