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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClassGit2
    {
        static int subdiscr1 = 4;// 4;// 6;
        static int discr1 = 2;// 3;//4;
                      // int discr2 dn xrhsimopoieitai
        static int subdiscr1_shell = 6;//14;
        static int discr1_shell = 1;
        static int graphene_sheets_number = 2; //periektikothta 0.525% 
        static double scale_factor = 1; //PROSOXH

        public static void /*(Model, double[])*/ RunExample()
        {
            int discr3 = discr1 * subdiscr1;// 23;
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = SeperateIntegrationClassCheck.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90;
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;
            
            bool run_new_corner = true; 
            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3DcornerGit(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);

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

            Dictionary<int, HashSet<INode>> cornerNodes = rveBuilder.cornerNodes;
            #endregion

            #region setup solver problem and initialize
            //Setup solver
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = 1E-5,
                MaxIterationsProvider = new FixedMaxIterationsProvider(100)
            };
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse()); //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();

            var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x => x.Item1, x => x.Value);

            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);
            var fetiSolverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
            fetiSolverBuilder.StiffnessDistribution =StiffnessDistributionType.HeterogeneousCondensed;
            fetiSolverBuilder.Preconditioning = new DirichletPreconditioning();
            FetiDPSolverSerial fetiSolver = fetiSolverBuilder.Build(model, cornerNodeSelection);
            model.ConnectDataStructures();
                     
            RunAnalysis(model, fetiSolver);
            #endregion



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


    }
}
