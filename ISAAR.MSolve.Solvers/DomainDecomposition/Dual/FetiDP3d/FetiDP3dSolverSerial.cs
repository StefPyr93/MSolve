using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MPI;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

//TODO: Add time logging
//TODO: Use a base class for the code that is identical between FETI-1 and FETI-DP.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d
{
    public class FetiDP3dSolverSerial : ISolverMpi, IFetiSolver
    {
        internal const string name = "FETI-DP Solver"; // for error messages and logging
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IFreeDofDisplacementsCalculator displacementsCalculator;
        private readonly DofOrderer dofOrderer;
        private readonly FetiDPDofSeparatorSerial dofSeparator;
        private readonly IFetiDPInterfaceProblemSolver interfaceProblemSolver;
        private readonly LagrangeMultipliersEnumeratorSerial lagrangesEnumerator;
        private readonly FetiDP3dMatrixManagerSerial matrixManager;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly IFetiPreconditionerFactory precondFactory;
        private readonly IFetiPreconditioningOperations preconditioning; //TODO: perhaps this should be hidden inside IFetiPreconditionerFactory
        private readonly IStiffnessDistribution stiffnessDistribution;
        private readonly FetiDPSubdomainGlobalMappingSerial subdomainGlobalMapping;

        private bool factorizeInPlace = true;
        private IFetiDPFlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private IFetiPreconditioner preconditioner;

        public Vector previousLambda { get; set; }

        public bool usePreviousLambda { get; set; }
        public FetiDP3dSolverSerial(IModel model, ICornerNodeSelection cornerNodeSelection, IMidsideNodesSelection midsideNodesSelection, IAugmentationConstraintsFactory augmentationConstraintsFactory,
            IFetiDP3dMatrixManagerFactory matrixManagerFactory, IFetiPreconditioningOperations preconditioning,
            ICrosspointStrategy crosspointStrategy, PcgSettings pcgSettings, StiffnessDistributionType stiffnessDistributionType,
            bool reorthogonalization)
        {
            this.msgHeader = $"{this.GetType().Name}: ";

            if (model.NumSubdomains == 1) throw new InvalidSolverException(msgHeader 
                + $"This solver cannot be used if there is only 1 subdomain");
            this.model = model;

            this.Logger = new SolverLoggerSerial(name);

            // Connectivity
            this.CornerNodes = cornerNodeSelection;
            this.dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            this.dofSeparator = new FetiDPDofSeparatorSerial(model, cornerNodeSelection);
            this.lagrangesEnumerator = new LagrangeMultipliersEnumeratorSerial(model, crosspointStrategy, dofSeparator);
            this.augmentationConstraints = augmentationConstraintsFactory.CreateAugmentationConstraints(model, 
                midsideNodesSelection, dofSeparator, lagrangesEnumerator);

            // Matrix managers and linear systems
            this.matrixManager = new FetiDP3dMatrixManagerSerial(model, this.dofSeparator, lagrangesEnumerator, 
                augmentationConstraints ,matrixManagerFactory);
            //TODO: This will call HandleMatrixWillBeSet() once for each subdomain. For now I will clear the data when 
            //      BuildMatrices() is called. Redesign this.
            //matrixManager.LinearSystem.MatrixObservers.Add(this); 

            // Preconditioning
            this.preconditioning = preconditioning;
            this.precondFactory = new FetiPreconditionerSerial.Factory();

            // Interface problem
            if (reorthogonalization)
            {
                this.interfaceProblemSolver = new FetiDP3dInterfaceProblemSolverReorthogonalizationSerial(
                    model, pcgSettings, augmentationConstraints);
            }
            else
            {
                this.interfaceProblemSolver = new FetiDP3dInterfaceProblemSolverSerial(
                    model, pcgSettings, augmentationConstraints);
            }
            this.displacementsCalculator = new FetiDP3dFreeDofDisplacementsCalculatorSerial(
                model, dofSeparator, lagrangesEnumerator,augmentationConstraints, matrixManager);

            // Homogeneous/heterogeneous problems
            // Homogeneous/heterogeneous problems
            if (stiffnessDistributionType == StiffnessDistributionType.Homogeneous)
            {
                this.stiffnessDistribution = new HomogeneousStiffnessDistributionSerial(model, dofSeparator,
                    new FetiDPHomogeneousDistributionLoadScaling(dofSeparator));
            }
            else if (stiffnessDistributionType == StiffnessDistributionType.HeterogeneousLumped)
            {
                this.stiffnessDistribution = new HeterogeneousLumpedStiffnessDistributionSerial(model, dofSeparator,
                    lagrangesEnumerator, matrixManager, new FetiDPHeterogeneousDistributionLoadScaling(dofSeparator));
            }
            else if (stiffnessDistributionType == StiffnessDistributionType.HeterogeneousCondensed)
            {
                this.stiffnessDistribution = new HeterogeneouspCondensedStiffnessDistributionSerial(model, dofSeparator,
                    lagrangesEnumerator, matrixManager, new FetiDPHeterogeneousDistributionLoadScaling(dofSeparator));
            }

            this.subdomainGlobalMapping = new FetiDPSubdomainGlobalMappingSerial(model, dofSeparator, stiffnessDistribution);
        }

        public ICornerNodeSelection CornerNodes { get; }
        public ISolverLogger Logger { get; }
        public string Name => name;
        public INodalLoadDistributor NodalLoadDistributor => stiffnessDistribution;

        public IFetiDPInterfaceProblemSolver InterfaceProblemSolver => interfaceProblemSolver;

        /// <summary>
        ///  builds Kff of each subdomain
        /// </summary>
        /// <param name="elementMatrixProvider"></param>
        public void BuildGlobalMatrix(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            Logger.StartMeasuringTime();

            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader
                        + $" Assembling the free-free stiffness matrix of subdomain {subdomain.ID}");
                    IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
                    subdomainMatrices.BuildFreeDofsMatrix(subdomain.FreeDofOrdering, elementMatrixProvider); 
                }
            }

            Logger.LogCurrentTaskDuration("Matrix assembly");

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
        }

        public Vector GatherGlobalDisplacements()
        {
            return subdomainGlobalMapping.GatherGlobalDisplacements(
                sub => matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.SolutionConcrete);
        }

        public ILinearSystemMpi GetLinearSystem(ISubdomain subdomain)
            => matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).LinearSystem;
        
        public void HandleMatrixWillBeSet()
        {
            isStiffnessModified = true;
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader + $"Clearing saved matrices of subdomain {subdomain.ID}.");
                    matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ClearMatrices();
                }
            }

            flexibility = null;
            preconditioner = null;
            matrixManager.ClearInverseCoarseProblemMatrix();
        }

        public void Initialize()
        {
            Logger.StartMeasuringTime();

            CornerNodes.Update();

            // Define the various dof groups
            dofSeparator.SeparateDofs(matrixManager);
            //FetiDPDofSeparationLogging.PrintDofSeparationSerial(model, dofSeparator);

            //TODO: B matrices could also be reused in some cases
            // Define lagrange multipliers and boolean matrices. 
            lagrangesEnumerator.CalcBooleanMatrices(dofSeparator.GetRemainderDofOrdering);

            //Define augmentation constraints and boolean matrices.
            augmentationConstraints.CalcAugmentationMappingMatrices();


            // Log dof statistics
            Logger.LogCurrentTaskDuration("Dof ordering");
            Logger.LogNumDofs("Lagrange multipliers", lagrangesEnumerator.NumLagrangeMultipliers);
            Logger.LogNumDofs("Corner dofs", dofSeparator.NumGlobalCornerDofs);
            Logger.LogNumDofs("AugmentationConstraints", augmentationConstraints.NumGlobalAugmentationConstraints);

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            stiffnessDistribution.Update();
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
        {
            Logger.StartMeasuringTime();

            // Order dofs
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.ConnectivityModified)
                {
                    matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).HandleDofOrderingWillBeModified(); //TODO: Not sure about this
                }
            }

            // This should not create subdomain-global mappings which require MPI communication
            //TODO: What about subdomain-global mappings, especially for boundary dofs? Who should create them? 
            dofOrderer.OrderFreeDofs(model);

            if (alsoOrderConstrainedDofs)
            {
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);
                }
            }

            // Log dof statistics
            Logger.LogCurrentTaskDuration("Dof ordering");
            Logger.LogNumDofs("Global dofs", model.GlobalDofOrdering.NumGlobalFreeDofs);
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            #region output
            if (CnstValues.printNRstiffnessMatrices && CnstValues.analyzerInfoIsSolutionForNRiters)
            {
                string print_path_gen = (new CnstValues()).exampleOutputPathGen + @"\subdomain_matrices_and_data\Subdomain{0}Iter{1}Loads.txt";
                string print_path_gen_st = (new CnstValues()).exampleOutputPathGen + @"\subdomain_matrices_and_data\Subdomain{0}Iter{1}Stiffness.txt";
                foreach (var subd in model.EnumerateSubdomains())
                {
                    var linearSystem = GetLinearSystem(subd); int subdId = subd.ID;

                    MatlabWriter eriter = new MatlabWriter();
                    //string print_path_gen = (new CnstValues()).exampleOutputPathGen + @"\subdomain_matrices_and_data\GlobalSuiteMatStress{0}LoadStep{1}Iter{2}_.txt";
                    //string print_path = string.Format(print_path_gen, CnstValues.stressIncrNo, CnstValues.analyzerLoadingStep, CnstValues.analyzerNRIter);

                    //var mat2 = DokColMajor.CreateFromDense(linearSystem.Matrix, 1e-14);

                    string print_path = string.Format(print_path_gen, subdId, CnstValues.analyzerNRIter);
                    eriter.WriteToFile(linearSystem.RhsVector, print_path, false);

                    bool print_stiffnesses = true;
                    if (print_stiffnesses)
                        {
                        string print_path_sti = string.Format(print_path_gen_st, subdId, CnstValues.analyzerNRIter);
                        eriter.WriteToFile((ISparseMatrix)GetLinearSystem(subd).Matrix, print_path_sti, false);
                    }
                    //string print_path_gen2 = (new CnstValues()).exampleOutputPathGen + @"\subdomain_matrices_and_data\GlobalSuiteRHSStress{0}LoadStep{1}Iter{2}_.txt";
                    //string print_path2 = string.Format(print_path_gen2, CnstValues.analyzerLoadingStep, CnstValues.analyzerNRIter);
                    //(new ISAAR.MSolve.LinearAlgebra.Output.Array1DWriter()).WriteToFile(linearSystem.RhsConcrete.CopyToArray(), print_path2);
                }
                
            }
            #endregion
            if (isStiffnessModified)
            {
                // Separate the stiffness matrix
                Logger.StartMeasuringTime();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ExtractCornerRemainderSubmatrices();
                    }
                }
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                // Calculate the preconditioner before factorizing each subdomain's Krr.
                // The inter-subdomain stiffness distribution may have changed even if a subdomain's stiffness is the same.
                Logger.StartMeasuringTime();
                if (preconditioning.ReorderInternalDofsForFactorization) dofSeparator.ReorderInternalDofs(matrixManager);
                preconditioner = precondFactory.CreatePreconditioner(preconditioning, model, dofSeparator, lagrangesEnumerator,
                    matrixManager, stiffnessDistribution);
                Logger.LogCurrentTaskDuration("Calculating preconditioner");



                // Factorize each subdomain's Krr
                Logger.StartMeasuringTime();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        //TODO: If I can reuse Krr, I can also reuse its factorization. Therefore this must be inPlace. In contrast, FETI-1 needs Kff intact for Stiffness distribution, in the current design).
                        Debug.WriteLine(msgHeader
                            + $"Inverting the remainder-remainder stiffness matrix of subdomain {subdomain.ID} in place.");
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).InvertKrr(true);
                    }
                }

                // Calculate FETI-DP coarse problem matrix
                matrixManager.CalcInverseCoarseProblemMatrix(CornerNodes);
                flexibility = new FetiDP3dFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, 
                    augmentationConstraints, matrixManager);
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                isStiffnessModified = false;
            }

            // Calculate FETI-DP coarse problem rhs 
            //TODO: rename this it is not coarse problem rhs only
            Logger.StartMeasuringTime();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ExtractCornerRemainderRhsSubvectors();
            }
            matrixManager.CalcCoarseProblemRhs();
            Logger.LogCurrentTaskDuration("Calculating coarse problem rhs");

            Logger.StartMeasuringTime();
            // Calculate the norm of the forces vector |f| = |K*u|. It is needed to check the convergence of PCG.
            double globalForcesNorm = globalForcesNorm = subdomainGlobalMapping.CalcGlobalForcesNorm(
                    sub => matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.RhsConcrete);

            if (CnstValues.printSolver_run_overwrite_data_region)
            {
                PrintLagrangeEqsData();
                PrintCornerEqsData();
            }

            // Solve interface problem
            InterfaceProblemSolver.PreviousLambda = previousLambda;
            InterfaceProblemSolver.UsePreviousLambda = usePreviousLambda;
            Vector lagranges = InterfaceProblemSolver.SolveInterfaceProblem(matrixManager, lagrangesEnumerator,
                flexibility, preconditioner, globalForcesNorm, Logger);
            if (usePreviousLambda) { previousLambda = lagranges; }
            
            Logger.LogCurrentTaskDuration("Solving interface problem");

            // Calculate the displacements of each subdomain
            Logger.StartMeasuringTime();
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);
            Vector globalU = GatherGlobalDisplacements();// sudomainDisplacements);
            ScatterGLobalSolutionU(matrixManager, model, globalU);
            Logger.LogCurrentTaskDuration("Calculate displacements from lagrange multipliers");

            Logger.IncrementAnalysisStep();
        }

        private void ScatterGLobalSolutionU(FetiDP3dMatrixManagerSerial matrixManager, IModel model, Vector globalU)
        {
            foreach (var entry in model.GlobalDofOrdering.GlobalFreeDofs)
            {
                var node = entry.row;
                var dofType = entry.col;
                foreach (var subdomaain in node.SubdomainsDictionary.Values)
                {
                    matrixManager.GetSubdomainMatrixManager(model.GetSubdomain(subdomaain.ID)).LinearSystem.SolutionConcrete[subdomaain.FreeDofOrdering.FreeDofs[node, dofType]] = globalU[entry.val];

                }
            }

        }

        private void PrintLagrangeEqsData()
        {
            List<int>[] subdomainLocals = new List<int>[model.EnumerateSubdomains().Count()];
            List<int>[] subdomainGobEqs = new List<int>[model.EnumerateSubdomains().Count()];
            List<int>[] subdomainPlMinus = new List<int>[model.EnumerateSubdomains().Count()];

            for (int i1 = 0; i1 < model.EnumerateSubdomains().Count(); i1++)
            {
                subdomainLocals[i1] = new List<int>();
                subdomainGobEqs[i1] = new List<int>();
                subdomainPlMinus[i1] = new List<int>();
            }

            for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; i++)
            {
                var lagrange = lagrangesEnumerator.LagrangeMultipliers[i];
                int subdId_pl = lagrange.SubdomainPlus.ID;
                int subdId_mn = lagrange.SubdomainMinus.ID;
                int dofId_local_plus = lagrange.SubdomainPlus.FreeDofOrdering.FreeDofs[lagrange.Node, lagrange.DofType];
                int dofId_local_minus = lagrange.SubdomainMinus.FreeDofOrdering.FreeDofs[lagrange.Node, lagrange.DofType];

                subdomainLocals[subdId_pl].Add(dofId_local_plus + 1);
                subdomainGobEqs[subdId_pl].Add(i + 1);
                subdomainPlMinus[subdId_pl].Add(+1);

                subdomainLocals[subdId_mn].Add(dofId_local_minus + 1);
                subdomainGobEqs[subdId_mn].Add(i + 1);
                subdomainPlMinus[subdId_mn].Add(-1);
            }

            var LagrangeNodeIds = new List<int>();
            var LagrangeNodeCouplingWaysNums = new Dictionary<int, int>();
            int previousNodeId = -1;
            int position = 0;
            for (int i = 0; i < lagrangesEnumerator.NumLagrangeMultipliers; i++)
            {
                int currentNodeId = lagrangesEnumerator.LagrangeMultipliers[i].Node.ID;
                if (!(currentNodeId == previousNodeId))
                {
                    LagrangeNodeIds.Add(currentNodeId);
                    position = LagrangeNodeIds.Count();
                    LagrangeNodeCouplingWaysNums[position] = 1;
                    previousNodeId = currentNodeId;


                }
                else
                {
                    LagrangeNodeCouplingWaysNums[position] = LagrangeNodeCouplingWaysNums[position] + 1;
                    //coupling_ways_num += 1;
                }
            }



            //string subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug\RVE_database\rve_no_1";
            //string subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_1";
            string subdomainOutputPath_gen = (new CnstValues()).exampleOutputPathGen+ @"\model_overwrite\subdomain_data_solver";


            for (int i1 = 0; i1 < model.EnumerateSubdomains().Count(); i1++)
            {
                var path1 = subdomainOutputPath_gen + @"\subdomainLagrangesLocals" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainLocals[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainLagrangesGlobalEqs" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainGobEqs[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainLagrangesPlMinus" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainPlMinus[i1], path1);

                path1 = subdomainOutputPath_gen + @"\RB_Nodes_IDs_MSOLVE_wise" + ".txt";
                WriteToFileVector(LagrangeNodeIds, path1);

                path1 = subdomainOutputPath_gen + @"\RB_Nodes_IDs_couplingwaysNum_MSOLVE_wise" + ".txt";
                WriteToFileVector(LagrangeNodeCouplingWaysNums.Values.ToList(), path1);
            }


        }

        private void PrintCornerEqsData()
        {
            List<int>[] subdomainLocals = new List<int>[model.EnumerateSubdomains().Count()];
            List<int>[] subdomainGobEqs = new List<int>[model.EnumerateSubdomains().Count()];

            for (int i1 = 0; i1 < model.EnumerateSubdomains().Count(); i1++)
            {
                subdomainLocals[i1] = new List<int>();
                subdomainGobEqs[i1] = new List<int>();

                var subdomain = model.EnumerateSubdomains().ElementAt(i1);

                foreach ((INode node, IDofType dof, _) in dofSeparator.GetCornerDofOrdering(subdomain) )
                {
                    subdomainLocals[i1].Add(subdomain.FreeDofOrdering.FreeDofs[node, dof] + 1);
                    subdomainGobEqs[i1].Add(dofSeparator.GlobalCornerDofOrdering[node, dof] + 1);
                }
            }

            //string subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug\RVE_database\rve_no_1";
            //string subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_1";
            string subdomainOutputPath_gen = (new CnstValues()).exampleOutputPathGen + @"\model_overwrite\subdomain_data_solver";


            for (int i1 = 0; i1 < model.EnumerateSubdomains().Count(); i1++)
            {
                var path1 = subdomainOutputPath_gen + @"\subdomainCornerLocals" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainLocals[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainCornerGlobalEqs" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainGobEqs[i1], path1);

            }

        }

        public static void WriteToFileVector(List<int> array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.Count(); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

            writer2.Dispose();

        }

        public class Builder
        {
            private readonly IFetiDP3dMatrixManagerFactory matrixManagerFactory;

            public Builder(IFetiDP3dMatrixManagerFactory matrixManagerFactory)
            {
                this.matrixManagerFactory = matrixManagerFactory;
            }

            public IAugmentationConstraintsFactory AugmentationConstraintsFactory { get; set; } = new AugmentationConstraints.Factory();
            
            public ICrosspointStrategy CrosspointStrategy { get; set; } = new FullyRedundantConstraints();

            public IDofOrderer DofOrderer { get; set; } =
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public PcgSettings PcgSettings { get; set; } = new PcgSettings();

            public IFetiPreconditioningOperations Preconditioning { get; set; } = new DirichletPreconditioning();

            public bool Reorthogonalization { get; set; } = false;

            public StiffnessDistributionType StiffnessDistribution { get; set; } = StiffnessDistributionType.Homogeneous;

            public FetiDP3dSolverSerial Build(IModel model, ICornerNodeSelection cornerNodeSelection,
                IMidsideNodesSelection midsideNodesSelection)
            {
                return new FetiDP3dSolverSerial(model, cornerNodeSelection, midsideNodesSelection, 
                    AugmentationConstraintsFactory, matrixManagerFactory, Preconditioning, CrosspointStrategy,
                    PcgSettings, StiffnessDistribution, Reorthogonalization);
            }
        }
    }
}
