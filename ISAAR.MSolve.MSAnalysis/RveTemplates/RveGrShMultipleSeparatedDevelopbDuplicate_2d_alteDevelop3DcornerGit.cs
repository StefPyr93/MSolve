using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.SamplesConsole.SupportiveClasses;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.MSAnalysis.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Creates a elastic matrix rve with embedded graphene sheets, with the asumption of damage behaviour for the material interface.
    /// Use of model separation methods is made.
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGit : IRVEbuilder //IdegenerateRVEbuilder
    {
        //origin: RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3D opou ginetai debug h 3D FETI me Ser
        //changes:  epiprostheta merika corner nodes


        public int[] hexaPrint { get; private set; }
        public int[] cohePrint { get; private set; }
        public int[] shellPrint { get; private set; }
        public Dictionary<int, double[]> CornerNodesIds { get; private set; }
        public Dictionary<int, int[]> CornerNodesIdAndsubdomains { get; private set; }
        public Dictionary<int, int[]> subdFreeBRNodes { get; private set; }
        public string subdomainOutputPath { get; private set; }
        public IList<Node> EmbeddedNodes { get; private set; }
        public Dictionary<ISubdomain, List<Node>> RveMatrixSubdomainInnerNodes { get; private set; }
        public List<List<int>> extraConstraintsNoeds { get; set; }
        public List<List<int>> extraConstraintsNoedsAve { get; set; }

        private bool decomposeModel;
        
        public Dictionary<int, HashSet<INode>> cornerNodes;
        public ISolver GetAppropriateSolver(Model model)
        {
            if (decomposeModel)
            {
                //Setup solver
                var pcgSettings = new PcgSettings()
                {
                    ConvergenceTolerance = 1E-5,
                    MaxIterationsProvider = new FixedMaxIterationsProvider(100)
                };
                //var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
                //interfaceSolverBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1);
                //interfaceSolverBuilder.PcgConvergenceTolerance = 1E-10;
                var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
                //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();

                var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x=>x.Item1,x=>x.Value);

                var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);
                var fetiSolverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
                //fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
                fetiSolverBuilder.StiffnessDistribution = StiffnessDistributionType.HeterogeneousCondensed;
                fetiSolverBuilder.Preconditioning = new DirichletPreconditioning();
                FetiDPSolverSerial fetiSolver = fetiSolverBuilder.Build(model,cornerNodeSelection);
                return (ISolver)fetiSolver;
            }
            else
            {
                return (new SkylineSolver.Builder()).BuildSolver(model);
            }
        }

        public Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
        int subdiscr1;
        int discr1;
        // int discr2 dn xrhsimopoieitai
        int discr3;
        int subdiscr1_shell;
        int discr1_shell;
        int graphene_sheets_number;
        public rveMatrixParameters mp;
        public grapheneSheetParameters gp;
        public string renumbering_vector_path;
        int RVE_id;
        public int[][] CornerNodesData;

        public RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGit(int RVE_id, bool decomposeModel, Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp,
            int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell, int graphene_sheets_number)
        {
            this.RVE_id = RVE_id;
            this.decomposeModel = decomposeModel;
            this.mpgp = mpgp;
            this.subdiscr1 = subdiscr1;
            this.discr1 = discr1;
            // int discr2 dn xrhsimopoieitai
            this.discr3 = discr3;
            this.subdiscr1_shell = subdiscr1_shell;
            this.discr1_shell = discr1_shell;
            this.graphene_sheets_number = graphene_sheets_number;

        }

        public IRVEbuilder Clone(int a) => new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGit(a, decomposeModel, mpgp, subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);

        public Tuple<Model, Dictionary<int, Node>, double> GetModelAndBoundaryNodes()
        {
            return Reference2RVEExample10000withRenumberingwithInput_forMS();
        }

        private Tuple<Model, Dictionary<int, Node>, double> Reference2RVEExample10000withRenumberingwithInput_forMS()
        {
            bool useInput = CnstValues.useInput_forRVE;
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain(1));
            Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();

            #region initial input format paths
            double[,] Dq;

            var rve_id_data = RVE_id.ToString();
            int path = 3;

            //renumbering_vector_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\REF_new_total_numbering.txt";
            if (path == 1)
            {
                renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
                renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);
            }
            if (path == 2)
            {
                renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2d_alte\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
                renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);
            }


            //string Fxk_p_komvoi_rve_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\Fxk_p_komvoi_rve.txt";
            string Fxk_p_komvoi_rve_path = null;
            if (path == 1)
            {
                Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\Fxk_p_komvoi_rve.txt";
                Fxk_p_komvoi_rve_path = string.Format(Fxk_p_komvoi_rve_path, rve_id_data);
            }
            if (path == 2)
            {
                Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2d_alte\RVE_database\rve_no_{0}\Fxk_p_komvoi_rve.txt";
                Fxk_p_komvoi_rve_path = string.Format(Fxk_p_komvoi_rve_path, rve_id_data);
            }



            //string o_xsunol_input_path_gen = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\o_xsunol_gs_";
            string o_xsunol_input_path_gen = null;
            if (path == 1) o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\o_xsunol_gs_";
            if (path == 2) o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2d_alte\RVE_database\rve_no_{0}\o_xsunol_gs_";
            if (path == 1 | path == 2)
            {
                o_xsunol_input_path_gen = string.Format(o_xsunol_input_path_gen, rve_id_data);
                o_xsunol_input_path_gen = o_xsunol_input_path_gen + "{0}.txt";
            }
            string subdomainOutputPath_gen = null;
            //if (path == 1) subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}";
            //if (path == 2) subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2d_alte\RVE_database\rve_no_{0}";
            var run_new_corner = false;
            if (run_new_corner)
            { //subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example3\input_matlab";
                subdomainOutputPath_gen = (new CnstValues()).exampleOutputPathGen ;
            }
            //{ subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_{0}"; }
            else
            { //subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example3\input_matlab";
                subdomainOutputPath_gen = (new CnstValues()).exampleOutputPathGen;
            }
            //{ subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_{0}"; }
            subdomainOutputPath = subdomainOutputPath_gen;
            if (path == 3) renumbering_vector_path = subdomainOutputPath_gen + @"\geometry_model\generated_ox_sunol_vecs\sunol_nodes_numbering.txt";
            if (path == 3) o_xsunol_input_path_gen = subdomainOutputPath_gen + @"\geometry_model\generated_ox_sunol_vecs\o_xsunol_generated_gs_{0}.txt";
            //subdomainOutputPath = string.Format(subdomainOutputPath_gen, rve_id_data);
            #endregion

            #region model parameters
            //int subdiscr1 = 6;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 23;
            //int subdiscr1_shell = 14;
            //int discr1_shell = 1;
            //mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
            gp = mpgp.Item2;


            //int graphene_sheets_number = 10;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];//TODO delete this
            double[][] ekk_xyz = new double[graphene_sheets_number][]; for (int i1 = 0; i1 < ekk_xyz.Length; i1++) { ekk_xyz[i1] = new double[3] { 0, 0, 0 }; };
            #endregion

            #region corner nodes data choose
            int[][] CornerNodesData;
            if (run_new_corner)
            { CornerNodesData = GetCornerNodesDataMoreCornerNodes(discr1, subdiscr1); }
            else
            { CornerNodesData = GetCornerNodesData(discr1, subdiscr1); }
            #endregion

            #region model geometry
            //bool useInput = true;
            bool run_debug = false;
            double b1 = 10; double b2 = 10; double sigma_f = 0.2;
            if (run_debug)
            { b1 = 0; b2 = 0; }

            double[][] o_xsunol_vectors = null;
            int[] sunol_nodes_numbering = null;
            int[] kanonas_renumbering_2 = null;


            if (!useInput)
            {
                //create geometry
                (o_xsunol_vectors, sunol_nodes_numbering, kanonas_renumbering_2) = CreateRandomGrSh_o_xsunol_and_renumbering(b1, b2, sigma_f,
                mp, ekk_xyz, model_o_x_parameteroi, run_debug);
                //print geometry
                string modelOutputPath_gen;
                if (!run_debug) modelOutputPath_gen = subdomainOutputPath_gen + @"\geometry_model\generated_ox_sunol_vecs";
                else modelOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_random_geometry_Msolve\REF2_50_000_renu_new_multiple_algorithms_check_develop_copy_for_progr_random_direct_in_C\generated\";
                RveDataPrintMethods.WriteModelDataOutput(modelOutputPath_gen, subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, mp, gp, kanonas_renumbering_2,
                    o_xsunol_vectors, sunol_nodes_numbering);
            }


            #region assignrenumbering
            renumbering renumbering;
            if (useInput) { renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path)); }
            else { renumbering = new renumbering(sunol_nodes_numbering); }
            #endregion

            extraConstraintsNoeds = GetExtraConstraintNodes2renumberingRandomData(discr1, subdiscr1, renumbering);
            extraConstraintsNoedsAve = GetExtraConstraintNodes2renumberingRandomDataOriginal(discr1, subdiscr1, renumbering);


            //TODO delete unesessary double arrays (Dq)
            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];


            if (useInput) { FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes); }
            else { FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering, boundaryNodes); }


            //domain separation ds1
            int totalSubdomains = discr1 * discr1 * discr1;
            if (decomposeModel) DdmCalculationsGeneral.BuildModelInterconnectionData(model);
            var decomposer = new AutomaticDomainDecomposer2(model, totalSubdomains);
            if (decomposeModel) decomposer.UpdateModel();
            var subdHexaIds = DdmCalculationsGeneral.DetermineHexaElementsSubdomainsFromModel(model);
            RveMatrixSubdomainInnerNodes = DdmCalculationsGeneral.DetermineRveSubdomainsInnerNodesFromModel(model);

            double volume = mp.L01 * mp.L02 * mp.L03;
            int hexaElementsNumber = model.ElementsDictionary.Count();

            //IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            //ds2
            int[] lowerCohesiveBound = new int[graphene_sheets_number]; int[] upperCohesiveBound = new int[graphene_sheets_number]; int[] grShElementssnumber = new int[graphene_sheets_number];
            for (int j = 0; j < graphene_sheets_number; j++)
            {
                if (useInput)
                {
                    string file_no = (j + 1).ToString();
                    string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                    FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_withRenumberingBondSlip(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                }
                else
                {
                    FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_from_MSOLVE_withRenumbering_from_MSOLVE(model, gp, renumbering, o_xsunol_vectors[j]);
                }
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                lowerCohesiveBound[j] = shellElementsNumber + element_counter_after_Adding_sheet + 1; //ds3
                upperCohesiveBound[j] = 2 * shellElementsNumber + element_counter_after_Adding_sheet;
                grShElementssnumber[j] = shellElementsNumber;
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }



            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = (EmbElementsIds).Select(x => model.ElementsDictionary[x]);
            //IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh


            //var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
            //var CohesiveGroupings = new EmbeddedCohesiveGrouping[EmbElementsIds.GetLength(0)];

            var hostSubGroups = new Dictionary<int, IEnumerable<Element>>();
            for (int i1 = 0; i1 < EmbElementsIds.GetLength(0); i1++)
            {
                if (useInput) { hostSubGroups.Add(EmbElementsIds[i1], FEMMeshBuilder.GetHostGroupForCohesiveElement(model.ElementsDictionary[EmbElementsIds[i1]], mp, model, renumbering_vector_path)); }
                else { hostSubGroups.Add(EmbElementsIds[i1], FEMMeshBuilder.GetHostGroupForCohesiveElement(model.ElementsDictionary[EmbElementsIds[i1]], mp, model)); }
                //hostSubGroups.Add(EmbElementsIds[i1], FEMMeshBuilder.GetHostGroupForCohesiveElement(model.ElementsDictionary[EmbElementsIds[i1]], mp, model, renumbering_vector_path));

                //var embeddedGroup_i1 = new List<Element>(1) { model.ElementsDictionary[EmbElementsIds[i1]] };
                //CohesiveGroupings[i1] = new EmbeddedCohesiveGrouping(model, hostGroup_i1, embeddedGroup_i1);
            }

            var CohesiveGroupping = new EmbeddedCohesiveSubGrouping(model, hostSubGroups, embdeddedGroup);

            #region corner node data generation


            if (!decomposeModel)
            {
                if (useInput) { (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model); }
                else { (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model, renumbering); }
            }


            #endregion


            //ds4
            if (decomposeModel)
            {
                //int[][] subdCohElementIds = DdmCalculationsGeneral.DetermineCoheiveELementsSubdomainsSimple(model, totalSubdomains);
                (int[][] subdCohElementIds, Dictionary<int, List<int>> reassignedHexas, Dictionary<int, int> hexaOriginalSubdomains,
                    Dictionary<int, List<int>> SubdomainNeedsHexas) =
                    DdmCalculationsGeneral.DetermineCoheiveELementsSubdomainsSimple_Alte4(model, totalSubdomains, lowerCohesiveBound,
                    upperCohesiveBound, grShElementssnumber);

                (Dictionary<int, List<int>> hexaAndTheirSharingSubdomains, Dictionary<int, int> DuplicateHexaOriginalSubdomain) =
                   DdmCalculationsAlterna2.GetHexaSharingSubdomains(SubdomainNeedsHexas, subdHexaIds, false);
                int[][] SubdomainNeedsHexasIds = DdmCalculationsPartb.ConvertIntListToArray(SubdomainNeedsHexas, totalSubdomains);

                (int[][] subdAdditionalHexaIds, Dictionary<int, int> NewHexaIdsAndTheirOriginals) = DdmCalculationsAlterna2.ModelAddDuplicateHexas(model, mp, hexaAndTheirSharingSubdomains, totalSubdomains);

                //int[][] subdHexaIdsNew = DdmCalculationsGeneral.ReassignHexas(subdHexaIds, reassignedHexas, hexaOriginalSubdomains);
                int[][] subdShellElementIds = DdmCalculationsGeneral.DetermineShellELementsSubdomains(model, totalSubdomains, subdCohElementIds,
                lowerCohesiveBound, upperCohesiveBound, grShElementssnumber);
                int[][] subdElementIds1 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdHexaIds, subdCohElementIds);
                int[][] subdElementIds2 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdElementIds1, subdShellElementIds);
                int[][] subdElementIds3 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdElementIds2, subdAdditionalHexaIds);
                DdmCalculationsGeneral.UndoModelInterconnectionDataBuild(model);
                DdmCalculations.SeparateSubdomains(model, subdElementIds3);

                model.ConnectDataStructures();
                bool isTrue = DdmCalculationsGeneral.CheckSubdomainsEmbeddingHostNodes(model, RveMatrixSubdomainInnerNodes);


                //PerfmormNesessaryChecksSubdomains(model);

                #region print extra data 

                if (useInput) { (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model); }
                else { (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model, renumbering); }

                #region find embedded
                EmbeddedNodes = new List<Node>();
                foreach (Element element in model.EnumerateElements())
                {
                    if (element.ID == 572)
                    {
                        string breakpoint = "here";
                    }

                    if (element.ElementType is IEmbeddedElement)
                    {
                        var element_I = element.ElementType as IEmbeddedElement;
                        foreach (EmbeddedNode embeddedNode in element_I.EmbeddedNodes)
                        {
                            if (!EmbeddedNodes.Contains(embeddedNode.Node))
                            { EmbeddedNodes.Add(embeddedNode.Node); }
                        }
                    }
                }
                #endregion

                DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);

                subdFreeBRNodes = new Dictionary<int, int[]>(); //nodeID, subdomainIDs  
                foreach (Node node in model.NodesDictionary.Values)
                {
                    if ((node.SubdomainsDictionary.Keys.Count() > 1 && !EmbeddedNodes.Contains(node)) && (!CornerNodesIds.Keys.Contains(node.ID) && (node.Constraints.Count() == 0)))
                    {
                        int[] subdIDs = node.SubdomainsDictionary.Keys.ToArray();
                        subdFreeBRNodes.Add(node.ID, subdIDs);
                    }
                }

                #endregion

                DdmCalculationsGeneral.UndoModelInterconnectionDataBuild(model);


                bool print_subdomain_data = CnstValues.print_subdomain_data;
                if (print_subdomain_data)
                {
                    //DdmCalculationsGeneral.PrintSubdomainDataForPostPro(subdHexaIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro(subdHexaIds, SubdomainNeedsHexasIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdFreeBRNodes, subdomainOutputPath, @"\subdomain_matrices_and_data\subdomainBRNodesAndSubd.txt");
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndsubdomains, subdomainOutputPath, @"\subdomain_matrices_and_data\CornerNodesAndSubdIds.txt");
                }

                bool get_subdomain_data = true;
                if (get_subdomain_data)
                {
                    (hexaPrint, cohePrint, shellPrint) = DdmCalculationsGeneral.GetSubdomainDataForPostPro(subdHexaIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                }

            }

            return new Tuple<Model, Dictionary<int, Node>, double>(model, boundaryNodes, volume);

        }

        private int[][] GetCornerNodesData(int discr1, int subdiscr1)
        {
            int total_corner_nodes_num = (discr1 - 1) * (discr1 - 1) * (discr1 - 1);
            int[][] CornerNodesData = new int[total_corner_nodes_num][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            int thesi = 0;
            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {
                    for (int i3 = 0; i3 < discr1 - 1; i3++)
                    {
                        CornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, (subdiscr1) + i2 * subdiscr1, (subdiscr1) + i3 * subdiscr1 };
                        thesi++;
                    }
                }
            }

            return CornerNodesData;
        }

        private int[][] GetCornerNodesDataMoreCornerNodes(int discr1, int subdiscr1)
        {
            int total_corner_nodes_num = (discr1 - 1) * (discr1 - 1) * (discr1 - 1);
            int[][] CornerNodesData = new int[total_corner_nodes_num][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            int thesi = 0;
            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {
                    for (int i3 = 0; i3 < discr1 - 1; i3++)
                    {
                        CornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, (subdiscr1) + i2 * subdiscr1, (subdiscr1) + i3 * subdiscr1 };
                        thesi++;
                    }
                }
            }

            int peripheral_corner_nodes_num = (discr1 - 1) * (discr1 - 1) * 6;
            int[][] peripheralCornerNodesData = new int[peripheral_corner_nodes_num][];
            thesi = 0;
            for (int i2 = 0; i2 < discr1 - 1; i2++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    peripheralCornerNodesData[thesi] = new int[3] { +1, (subdiscr1) + i2 * subdiscr1, (subdiscr1) + i3 * subdiscr1 };
                    thesi++;
                }
            }

            for (int i2 = 0; i2 < discr1 - 1; i2++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    peripheralCornerNodesData[thesi] = new int[3] { discr1 * subdiscr1 - 1, (subdiscr1) + i2 * subdiscr1, (subdiscr1) + i3 * subdiscr1 };
                    thesi++;
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {

                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    peripheralCornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, +1, (subdiscr1) + i3 * subdiscr1 };
                    thesi++;
                }

            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {

                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    peripheralCornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, discr1 * subdiscr1 - 1, (subdiscr1) + i3 * subdiscr1 };
                    thesi++;
                }

            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {

                    peripheralCornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, (subdiscr1) + i2 * subdiscr1, +1 };
                    thesi++;

                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {

                    peripheralCornerNodesData[thesi] = new int[3] { (subdiscr1) + i1 * subdiscr1, (subdiscr1) + i2 * subdiscr1, discr1 * subdiscr1 - 1 };
                    thesi++;

                }
            }

            int[][] moreCornerNodeData = new int[20][];
            moreCornerNodeData[0] = new int[3] { 7, 4, 7 };
            moreCornerNodeData[1] = new int[3] { 1, 4, 7 };
            moreCornerNodeData[2] = new int[3] { 7, 4, 1 };
            moreCornerNodeData[3] = new int[3] { 1, 4, 1 };
            moreCornerNodeData[4] = new int[3] { 4, 1, 7 };
            moreCornerNodeData[5] = new int[3] { 4, 7, 7 };
            moreCornerNodeData[6] = new int[3] { 4, 1, 1 };
            moreCornerNodeData[7] = new int[3] { 4, 7, 1 };
            moreCornerNodeData[8] = new int[3] { 7, 1, 4 };
            moreCornerNodeData[9] = new int[3] { 7, 7, 4 };
            moreCornerNodeData[10] = new int[3] { 1, 1, 4 };
            moreCornerNodeData[11] = new int[3] { 1, 7, 4 };

            moreCornerNodeData[12] = new int[3] { 1, 1, 1 };
            moreCornerNodeData[13] = new int[3] { 1, 1, 7 };
            moreCornerNodeData[14] = new int[3] { 1, 7, 7 };
            moreCornerNodeData[15] = new int[3] { 1, 7, 1 };
            moreCornerNodeData[16] = new int[3] { 7, 1, 1 };
            moreCornerNodeData[17] = new int[3] { 7, 1, 7 };
            moreCornerNodeData[18] = new int[3] { 7, 7, 7 };
            moreCornerNodeData[19] = new int[3] { 7, 7, 1 };



            var cornersList1 = CornerNodesData.ToList();
            var cornersList2 = peripheralCornerNodesData.ToList();
            var moreCornerList = moreCornerNodeData.ToList();
            cornersList1.AddRange(cornersList2);
            cornersList1.AddRange(moreCornerList);

            return cornersList1.ToArray();
        }

        private void PerfmormNesessaryChecksSubdomains(Model model)
        {
            int subdToCheck = 41;

            List<Node> TotalNodes = new List<Node>();
            List<Node> TotalNodes2 = new List<Node>();
            List<Node> Hexa8Nodes = new List<Node>();
            List<Node> Hexa8Nodes2 = new List<Node>();

            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements.Values)
            {
                foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                {
                    TotalNodes.Add(node);
                }
            }
            TotalNodes = TotalNodes.Union(TotalNodes2).ToList();
            int ch01 = TotalNodes.Count();


            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements.Values)
            {
                if (element.ElementType is Hexa8NonLinear)
                {
                    foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                    {
                        Hexa8Nodes.Add(node);
                    }
                }
            }
            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements.Values)
            {
                if (element.ElementType is Shell8NonLinear)
                {
                    foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                    {
                        Hexa8Nodes2.Add(node);
                    }
                }
            }
            Hexa8Nodes = Hexa8Nodes.Union(Hexa8Nodes2).ToList();
            int ch020 = Hexa8Nodes.Count();

            var unconnected = TotalNodes.RemoveAll(x => Hexa8Nodes.Contains(x));

            var elements1 = TotalNodes[0].ElementsDictionary.Values.ToList();
            elements1.RemoveAll(x => (!(x.Subdomain.ID == 41)));
            elements1.RemoveAll(x => (!(x.ElementType is Hexa8NonLinear)));
            var elements2 = TotalNodes[1].ElementsDictionary.Values.ToList();
            elements2.RemoveAll(x => (!(x.Subdomain.ID == 41)));
            elements2.RemoveAll(x => (!(x.ElementType is Hexa8NonLinear)));

            var elemtns = model.ElementsDictionary.Values.ToList();
            elemtns.RemoveAll(x => (!CheckElementSameNodes(x, model.ElementsDictionary[4182])));

            var embElement2 = elemtns[2].ElementType as IEmbeddedElement;
            int[] subdIds2 = new int[embElement2.EmbeddedNodes.Count()];
            int[] nodesIDs2 = new int[embElement2.EmbeddedNodes.Count()];
            int[] elmntsIDs2 = new int[embElement2.EmbeddedNodes.Count()];
            for (int i1 = 0; i1 < embElement2.EmbeddedNodes.Count(); i1++)
            {
                subdIds2[i1] = model.ElementsDictionary[embElement2.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID].Subdomain.ID;
                //embElement2.EmbeddedNodes.ElementAt(0).EmbeddedInElement.Subdomain.ID dinei null apotelesma ara bug
                nodesIDs2[i1] = embElement2.EmbeddedNodes.ElementAt(i1).Node.ID;
                elmntsIDs2[i1] = embElement2.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID;
            }

            var embElement1 = elemtns[1].ElementType as IEmbeddedElement;
            int[] subdIds1 = new int[embElement1.EmbeddedNodes.Count()];
            int[] nodesIDs1 = new int[embElement1.EmbeddedNodes.Count()];
            int[] elmntsIDs1 = new int[embElement1.EmbeddedNodes.Count()];
            for (int i1 = 0; i1 < embElement1.EmbeddedNodes.Count(); i1++)
            {
                subdIds1[i1] = model.ElementsDictionary[embElement1.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID].Subdomain.ID;
                nodesIDs1[i1] = embElement1.EmbeddedNodes.ElementAt(i1).Node.ID;
                elmntsIDs1[i1] = embElement1.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID;
            }

            List<int> elmnt1AssemblyNodes = elemtns[1].ElementType.DofEnumerator.GetNodesForMatrixAssembly(elemtns[1]).ToList().ConvertAll(x => x.ID);
            List<int> elmnt2AssemblyNodes = elemtns[2].ElementType.DofEnumerator.GetNodesForMatrixAssembly(elemtns[2]).ToList().ConvertAll(x => x.ID);

            bool ch03 = elmnt1AssemblyNodes.Contains(2643);
            //false
            bool ch04 = elmnt1AssemblyNodes.Contains(2666);
            //false
            bool ch05 = elmnt2AssemblyNodes.Contains(2666);
            //true
            bool ch06 = elmnt2AssemblyNodes.Contains(2643);
            //true


        }

        public bool CheckElementSameNodes(Element element1, Element element2)
        {
            bool have8FirstSameNodes = true;
            for (int i1 = 0; i1 < 8; i1++)
            {
                if (!(element1.Nodes.ElementAt(i1).ID == element2.Nodes.ElementAt(i1).ID))
                {
                    have8FirstSameNodes = false;
                }
            }
            return have8FirstSameNodes;
        }

        private List<List<int>> GetExtraConstraintNodes2renumberingRandomDataOriginal(int discr1, int subdiscr1, renumbering renumbering)
        {
            List<List<int>> extraConstraintsNoeds = new List<List<int>>();
            //renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int i2 = 0; i2 < discr1 - 1; i2++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    int k2 = (subdiscr1) + i2 * subdiscr1;
                    int k3 = (subdiscr1) + i3 * subdiscr1;

                    for (int j1 = 0; j1 < discr1; j1++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        for (int i1 = 2 + j1 * subdiscr1; i1 < 1 + (1 + j1) * subdiscr1; i1++)
                        {
                            int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, k2 + 1, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                            oneEdgeNodes.Add(nodeID);
                        }


                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    int k1 = (subdiscr1) + i1 * subdiscr1;
                    int k3 = (subdiscr1) + i3 * subdiscr1;

                    for (int j2 = 0; j2 < discr1; j2++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        for (int i2 = 2 + j2 * subdiscr1; i2 < 1 + (1 + j2) * subdiscr1; i2++)
                        {
                            int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, i2, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                            oneEdgeNodes.Add(nodeID);
                        }
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {
                    int k1 = (subdiscr1) + i1 * subdiscr1;
                    int k2 = (subdiscr1) + i2 * subdiscr1;

                    for (int j3 = 0; j3 < discr1; j3++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        for (int i3 = 2 + j3 * subdiscr1; i3 < 1 + (1 + j3) * subdiscr1; i3++)
                        {
                            int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, k2 + 1, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                            oneEdgeNodes.Add(nodeID);
                        }
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            return extraConstraintsNoeds;
        }

        //midside node sta mesa twn plevrwn
        private List<List<int>> GetExtraConstraintNodes2renumberingRandomData(int discr1, int subdiscr1, renumbering renumbering)
        {
            List<List<int>> extraConstraintsNoeds = new List<List<int>>();
            //renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int i2 = 0; i2 < discr1 - 1; i2++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    int k2 = (subdiscr1) + i2 * subdiscr1;
                    int k3 = (subdiscr1) + i3 * subdiscr1;

                    for (int j1 = 0; j1 < discr1; j1++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j1 * subdiscr1;
                        int akri2 = 1 + (1 + j1) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i1 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, k2 + 1, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        oneEdgeNodes.Add(nodeID);
                        //for (int i1 = 2 + j1 * subdiscr1; i1 < 1 + (1 + j1) * subdiscr1; i1++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, k2 + 1, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}


                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    int k1 = (subdiscr1) + i1 * subdiscr1;
                    int k3 = (subdiscr1) + i3 * subdiscr1;

                    for (int j2 = 0; j2 < discr1; j2++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j2 * subdiscr1;
                        int akri2 = 1 + (1 + j2) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i2 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, i2, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        oneEdgeNodes.Add(nodeID);
                        //for (int i2 = 2 + j2 * subdiscr1; i2 < 1 + (1 + j2) * subdiscr1; i2++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, i2, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {
                    int k1 = (subdiscr1) + i1 * subdiscr1;
                    int k2 = (subdiscr1) + i2 * subdiscr1;

                    for (int j3 = 0; j3 < discr1; j3++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j3 * subdiscr1;
                        int akri2 = 1 + (1 + j3) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i3 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, k2 + 1, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        oneEdgeNodes.Add(nodeID);
                        //for (int i3 = 2 + j3 * subdiscr1; i3 < 1 + (1 + j3) * subdiscr1; i3++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, k2 + 1, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            return extraConstraintsNoeds;
        }
        
        //midside node sta mesa twn faces(mallon)
        private List<List<int>> GetExtraConstraintNodes2renumberingRandomDataV2(int discr1, int subdiscr1, renumbering renumbering)
        {
            List<List<int>> extraConstraintsNoeds = new List<List<int>>();
            //renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int j2 = 0; j2 < discr1; j2++)
            {
                //int k2 = (subdiscr1) + i2 * subdiscr1;
                int akri1_ = 1 + j2 * subdiscr1;
                int akri2_ = 1 + (1 + j2) * subdiscr1;
                int delta_ = akri2_ - akri1_;
                int miso_delta_ = (int)Math.Truncate((double)(delta_ / 2));

                int i2 = akri1_ + miso_delta_;


                for (int i3 = 0; i3 < discr1 - 1; i3++)
                {
                    int k3 = (subdiscr1) + i3 * subdiscr1;

                    for (int j1 = 0; j1 < discr1; j1++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j1 * subdiscr1;
                        int akri2 = 1 + (1 + j1) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i1 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, i2, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        oneEdgeNodes.Add(nodeID);
                        //for (int i1 = 2 + j1 * subdiscr1; i1 < 1 + (1 + j1) * subdiscr1; i1++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, k2 + 1, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}


                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int i1 = 0; i1 < discr1 - 1; i1++)
            {
                int k1 = (subdiscr1) + i1 * subdiscr1;

                for (int j3 = 0; j3 < discr1; j3++)
                {
                    //int k3 = (subdiscr1) + i3 * subdiscr1;
                    int akri1_ = 1 + j3 * subdiscr1;
                    int akri2_ = 1 + (1 + j3) * subdiscr1;
                    int delta_ = akri2_ - akri1_;
                    int miso_delta_ = (int)Math.Truncate((double)(delta_ / 2));

                    int i3 = akri1_ + miso_delta_;


                    for (int j2 = 0; j2 < discr1; j2++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j2 * subdiscr1;
                        int akri2 = 1 + (1 + j2) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i2 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, i2, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        oneEdgeNodes.Add(nodeID);
                        //for (int i2 = 2 + j2 * subdiscr1; i2 < 1 + (1 + j2) * subdiscr1; i2++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, i2, k3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            for (int j1 = 0; j1 < discr1; j1++)
            {
                //int k1 = (subdiscr1) + i1 * subdiscr1;
                int akri1_ = 1 + j1 * subdiscr1;
                int akri2_ = 1 + (1 + j1) * subdiscr1;
                int delta_ = akri2_ - akri1_;
                int miso_delta_ = (int)Math.Truncate((double)(delta_ / 2));

                int i1 = akri1_ + miso_delta_;

                for (int i2 = 0; i2 < discr1 - 1; i2++)
                {

                    int k2 = (subdiscr1) + i2 * subdiscr1;

                    for (int j3 = 0; j3 < discr1; j3++)
                    {
                        List<int> oneEdgeNodes = new List<int>();

                        int akri1 = 1 + j3 * subdiscr1;
                        int akri2 = 1 + (1 + j3) * subdiscr1;
                        int delta = akri2 - akri1;
                        int miso_delta = (int)Math.Truncate((double)(delta / 2));

                        int i3 = akri1 + miso_delta;
                        int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(i1, k2 + 1, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        oneEdgeNodes.Add(nodeID);
                        //for (int i3 = 2 + j3 * subdiscr1; i3 < 1 + (1 + j3) * subdiscr1; i3++)
                        //{
                        //    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(k1 + 1, k2 + 1, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        //    oneEdgeNodes.Add(nodeID);
                        //}
                        extraConstraintsNoeds.Add(oneEdgeNodes);
                    }
                }
            }

            return extraConstraintsNoeds;
        }

        // PROSOXH DEN ARKEI MONO TO PARAKATW NA GINEI UNCOMMENT WSTE NA GINEI IMPLEMENT TO IDegenerateRVEBuilder 
        //xreiazetai kai na xrhsimopoithei h katallhlh methodos tou femmeshbuilder gia to model and boundary nodes na dinei mono ta peripheral
        //public Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model model)
        //{
        //    return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
        //    //TODO:  Pithanws na epistrefetai apo GetModelAndBoundaryNodes ... AndConstraints.
        //}
        private void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        private Dictionary<int, HashSet<INode>> DefineCornerNodesPerSubdomainAndOtherwise(Dictionary<int, int[]> CornerNodesIdAndsubdomains, Model model)
        {
            // a copy of this is used in modelCreatorInput class
            Dictionary<int, HashSet<INode>> cornerNodesList = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());
            Dictionary<int, HashSet<INode>> cornerNodes = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());

            foreach (Subdomain subdomain in model.EnumerateSubdomains())
            {
                cornerNodesList.Add(subdomain.ID, new HashSet<INode>());
            }

            foreach (int CornerNodeID in CornerNodesIdAndsubdomains.Keys)
            {
                Node node1 = model.NodesDictionary[CornerNodeID];
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            //foreach (Subdomain subdomain in model.Subdomains)
            //{
            //    cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID]);
            //}

            return cornerNodesList;
        }

        private (Dictionary<int, double[]> CornerNodesIds, Dictionary<int, int[]> CornerNodesIdAndsubdomains, Dictionary<int, HashSet<INode>> cornerNodes)
            DefineCornerNodesFromCornerNodeData(int[][] CornerNodesData, Model model)
        {
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            CornerNodesIds = new Dictionary<int, double[]>(CornerNodesData.Length); //nodeID, coordinate data
            CornerNodesIdAndsubdomains = new Dictionary<int, int[]>(CornerNodesData.Length);//nodeID, subdIds            
                                                                                            //var CornerNodesSubdomains = new Dictionary<int, int[]> (CornerNodesData.Length); //nodeID, subdomains opou anhkei
            for (int i1 = 0; i1 < CornerNodesData.Length; i1++)
            {
                int h1 = CornerNodesData[i1][0]; int h2 = CornerNodesData[i1][1]; int h3 = CornerNodesData[i1][2];
                int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                double nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                double nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                double nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                double[] coordinates = new double[3] { nodeCoordX, nodeCoordY, nodeCoordZ };
                CornerNodesIds.Add(nodeID, coordinates);
                CornerNodesIdAndsubdomains.Add(nodeID, model.NodesDictionary[nodeID].SubdomainsDictionary.Keys.ToArray());
            }

            cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(CornerNodesIdAndsubdomains, model);

            return (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes);
        }

        private (Dictionary<int, double[]> CornerNodesIds, Dictionary<int, int[]> CornerNodesIdAndsubdomains, Dictionary<int, HashSet<INode>> cornerNodes)
            DefineCornerNodesFromCornerNodeData(int[][] CornerNodesData, Model model, renumbering renumbering)
        {
            //renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            CornerNodesIds = new Dictionary<int, double[]>(CornerNodesData.Length); //nodeID, coordinate data
            CornerNodesIdAndsubdomains = new Dictionary<int, int[]>(CornerNodesData.Length);//nodeID, subdIds            
                                                                                            //var CornerNodesSubdomains = new Dictionary<int, int[]> (CornerNodesData.Length); //nodeID, subdomains opou anhkei
            for (int i1 = 0; i1 < CornerNodesData.Length; i1++)
            {
                int h1 = CornerNodesData[i1][0]; int h2 = CornerNodesData[i1][1]; int h3 = CornerNodesData[i1][2];
                int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                double nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                double nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                double nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                double[] coordinates = new double[3] { nodeCoordX, nodeCoordY, nodeCoordZ };
                CornerNodesIds.Add(nodeID, coordinates);
                CornerNodesIdAndsubdomains.Add(nodeID, model.NodesDictionary[nodeID].SubdomainsDictionary.Keys.ToArray());
            }

            cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(CornerNodesIdAndsubdomains, model);

            return (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes);
        }

        public void UpdateStochasticCoefficientsProvider(IStochasticCoefficientsProvider2D coefficientsProvider)
        {
            coefficientsProvider.RandomVariables = null;
        }

        public (double[][] o_xsunol_vectors, int[] sunol_nodes_numbering, int[] kanonas_renumbering_2) CreateRandomGrSh_o_xsunol_and_renumbering(double b1, double b2, double sigma_f,
            rveMatrixParameters mp, double[][] ekk_xyz, o_x_parameters[] model_o_x_parameteroi, bool run_debug)
        {
            double[][] o_x_rve = FEMMeshBuilder.Build_o_x_rve_Coordinates(mp);
            IList<IStochasticCoefficientsProvider2D> coefficientsProviders =
                new List<IStochasticCoefficientsProvider2D> { new SpectralRepresentation2DRandomField(b1, b2, sigma_f, 0.01, 2 * Math.PI / ((double)20), 20) };
            for (int j = 0; j < graphene_sheets_number - 1; j++)
            { coefficientsProviders.Add(new SpectralRepresentation2DRandomField(b1, b2, sigma_f, 0.01, 2 * Math.PI / ((double)20), 20)); }

            double[][] o_xsunol_vectors = new double[graphene_sheets_number][];//mporei na xrhsimpopoiithei kai to ox_sunol_BUilder tou RveExamples builder tou pio prosfatou input commit
                                                                               //dld xwris to random 
            for (int j = 0; j < graphene_sheets_number; j++)
            {
                UpdateStochasticCoefficientsProvider(coefficientsProviders[j]); //origin: RandomGrapheneModelBuilder.FewElementsRVECheckExample2GrapheneSheets (Input 3 Ger Ody)
                int new_rows = 2 * gp.elem1 + 1;
                int new_lines = 2 * gp.elem2 + 1; //TODO taktopoihsh afta edw ta variables
                o_xsunol_vectors[j] = FEMMeshBuilder.ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, gp.L1, gp.L2, gp.elem1, gp.elem2, gp.a1_shell, ekk_xyz[j],
                model_o_x_parameteroi[j], coefficientsProviders[j]);// origin: RandomGrapheneModelBuilder.AddGrapheneSheet_with_o_x_parameters(......, IStochasticCoefficientsProvider2D coefficientsProvider)
            }

            //create random data for geom (origin:RVEkanoninkhsGewmetriasBuilder.Reference2RVEExample50_000withRenumberingwithInputFromMSOLVE() %637 apo input 4)
            sigma_f = 0.2; // apo to arxeio create_random_data_for_geom_programing_in_C tou fakelou tou parakatw rand data vec path
            string rand_data_vec_path = (new CnstValues()).rand_data_vec_path;
            savedRandomDataClass a = new savedRandomDataClass(PrintUtilities.ReadVector(rand_data_vec_path));
            Tuple<double[], double[], double[][]> RandomDataForGeomGiaSugkekrimenoRand;

            if (run_debug) { RandomDataForGeomGiaSugkekrimenoRand = RandomOrientations.CreateRandomDataForGeom(graphene_sheets_number, gp, mp, sigma_f, a); }
            else { RandomDataForGeomGiaSugkekrimenoRand = RandomOrientations.CreateRandomDataForGeom(graphene_sheets_number, gp, mp, sigma_f); }
            //return new Tuple<double[], double[], double[][]>(rot_phi_1, rot_phi_2, ekk_xyz);
            double[] rot_phi_1 = RandomDataForGeomGiaSugkekrimenoRand.Item1;
            double[] rot_phi_2 = RandomDataForGeomGiaSugkekrimenoRand.Item2;
            ekk_xyz = RandomDataForGeomGiaSugkekrimenoRand.Item3;

            //update o_xsunol_vectors  for rotation and translation
            for (int j = 0; j < graphene_sheets_number; j++)
            {
                o_xsunol_vectors[j] = RandomOrientations.modify_ox_sunol_forRotationAndTranslation(o_xsunol_vectors[j], rot_phi_1[j], rot_phi_2[j], ekk_xyz[j]);
            }
            #endregion

            #region create renumbering
            (int[] sunol_nodes_numbering, int[] kanonas_renumbering_2) = FEMMeshBuilder.GetTotalModelRenumbering(o_x_rve, o_xsunol_vectors, mp);
            #endregion

            return (o_xsunol_vectors, sunol_nodes_numbering, kanonas_renumbering_2);
        }

    }
}
