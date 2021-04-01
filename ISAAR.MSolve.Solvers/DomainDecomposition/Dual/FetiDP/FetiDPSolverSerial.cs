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
using System.Reflection;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;

//TODO: Add time logging
//TODO: Use a base class for the code that is identical between FETI-1 and FETI-DP.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPSolverSerial : ISolverMpi, IFetiSolver
    {
        internal const string name = "FETI-DP Solver"; // for error messages and logging
        private readonly IFreeDofDisplacementsCalculator displacementsCalculator;
        private readonly DofOrderer dofOrderer;
        private readonly FetiDPDofSeparatorSerial dofSeparator;
        private readonly FetiDPInterfaceProblemSolverSerial interfaceProblemSolver;
        private readonly LagrangeMultipliersEnumeratorSerial lagrangesEnumerator;
        private readonly FetiDPMatrixManagerSerial matrixManager;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly IFetiPreconditionerFactory precondFactory;
        private readonly IFetiPreconditioningOperations preconditioning; //TODO: perhaps this should be hidden inside IFetiPreconditionerFactory
        private readonly IStiffnessDistribution stiffnessDistribution;
        private readonly FetiDPSubdomainGlobalMappingSerial subdomainGlobalMapping;

        private bool factorizeInPlace = true;
        private FetiDPFlexibilityMatrixSerial flexibility;
        private bool isStiffnessModified = true;
        private IFetiPreconditioner preconditioner;

        public Vector previousLambda { get; set; } // TODO: implment lambda recycling for fetidpsolver as well
        public bool usePreviousLambda { get; set; }
        public FetiDPSolverSerial(IModel model, ICornerNodeSelection cornerNodeSelection,
            IFetiDPMatrixManagerFactory matrixManagerFactory, IFetiPreconditioningOperations preconditioning,
            ICrosspointStrategy crosspointStrategy, PcgSettings pcgSettings, StiffnessDistributionType stiffnessDistributionType)
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

            // Matrix managers and linear systems
            this.matrixManager = new FetiDPMatrixManagerSerial(model, this.dofSeparator, matrixManagerFactory);
            //TODO: This will call HandleMatrixWillBeSet() once for each subdomain. For now I will clear the data when 
            //      BuildMatrices() is called. Redesign this.
            //matrixManager.LinearSystem.MatrixObservers.Add(this); 

            // Preconditioning
            this.preconditioning = preconditioning;
            this.precondFactory = new FetiPreconditionerSerial.Factory();

            // Interface problem
            this.interfaceProblemSolver = new FetiDPInterfaceProblemSolverSerial(model, pcgSettings);
            this.displacementsCalculator = new FreeDofDisplacementsCalculatorSerial(model, dofSeparator, matrixManager,
                lagrangesEnumerator);

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

        public IFetiDPInterfaceProblemSolver InterfaceProblemSolver => interfaceProblemSolver;

        public ISolverLogger Logger { get; }
        public string Name => name;
        public INodalLoadDistributor NodalLoadDistributor => stiffnessDistribution;

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

            // Log dof statistics
            Logger.LogCurrentTaskDuration("Dof ordering");
            Logger.LogNumDofs("Lagrange multipliers", lagrangesEnumerator.NumLagrangeMultipliers);
            Logger.LogNumDofs("Corner dofs", dofSeparator.NumGlobalCornerDofs);

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
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ISingleSubdomainLinearSystemMpi linearSystem = 
                    matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).LinearSystem;
                linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            }

            if (isStiffnessModified)
            {
                // Separate the stiffness matrix 
                Logger.StartMeasuringTime();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ExtractCornerRemainderSubmatrices();
                        #region debug1
                        // TODOGer: edw ginetai extract Kcc, Krc kai Krr. en telei sto FetiDPsubdomainMatrixmanagerSkyline %108
                        #endregion
                    }
                }
                #region debug1
                /*
                var nodeCoords = new double[3] { -22.5, -22.5, -22.5 };
                int nodeId = ((Model)model).NodesDictionary.Values.Where(x => ((x.X1 == nodeCoords[0]) && (x.X2 == nodeCoords[1]) && (x.X3 == nodeCoords[2]))).ToList().ElementAt(0).ID;
                int subdID = 2; // ((Model)model).NodesDictionary[nodeId].SubdomainsDictionary.ElementAt(0).Key;

                var nodeCoords2 = new double[3][] { new double[3] { 0, -22.5, -22.5 }, new double[3] { -22.5,0, -22.5 }, new double[3] { -22.5, -22.5,0 } };
                int[] nodeIds  = nodeCoords2.Select(x=>((Model)model).NodesDictionary.Values.Where(y => ((y.X1 == x[0]) && y.X2 == x[1]) && (y.X3 == x[2])).ToList().ElementAt(0).ID).ToArray();

                var crossPointCoords = new double[3][] {  new double[3] { -22.5, 0, 0 }, new double[3] { 0, -22.5, 0 }, new double[3] { 0,0, -22.5 } };
                int[] crossPointIds = crossPointCoords.Select(x => ((Model)model).NodesDictionary.Values.Where(y => ((y.X1 == x[0]) && y.X2 == x[1]) && (y.X3 == x[2])).ToList().ElementAt(0).ID).ToArray();

                var cornerNodeCoords = new double[3] { 0, 0, 0 };
                int cornerNodeId = ((Model)model).NodesDictionary.Values.Where(x => ((x.X1 == cornerNodeCoords[0]) && (x.X2 == cornerNodeCoords[1]) && (x.X3 == cornerNodeCoords[2]))).ToList().ElementAt(0).ID;
                // Use reflection to set the necessary matrices
                var ch01 = matrixManager.GetFetiDPSubdomainMatrixManager(model.GetSubdomain(subdID));
                FieldInfo fi;
                fi = typeof(FetiDPSubdomainMatrixManagerSkyline).GetField("Krr", BindingFlags.NonPublic | BindingFlags.Instance);
                SkylineMatrix Krr = (SkylineMatrix)(fi.GetValue(ch01));

                fi = typeof(FetiDPSubdomainMatrixManagerSkyline).GetField("Krc", BindingFlags.NonPublic | BindingFlags.Instance);
                CscMatrix Krc = (CscMatrix)(fi.GetValue(ch01));

                //var dofSeparator = ((FetiDPFlexibilityMatrixSerial)flexibility).dofSeparator;
                DofTable subdRemainderDofs = dofSeparator.GetRemainderDofOrdering(model.GetSubdomain(subdID));
                DofTable subdCornerDofs = dofSeparator.GetCornerDofOrdering(model.GetSubdomain(subdID));

                int cornerNode_xdof_order = subdCornerDofs[model.GetNode(cornerNodeId), StructuralDof.TranslationX];
                //int node2_remainder_xdof_order = subdRemainderDofs[model.GetNode(nodeIds[0]), StructuralDof.TranslationX]; // TODO 1
                int crosspoint_remainder_xdof_order = subdRemainderDofs[model.GetNode(crossPointIds[0]), StructuralDof.TranslationX];
                             
                //double krc_value = Krc[node2_remainder_xdof_order, cornerNode_xdof_order];
                //double krr_value = Krr[node2_remainder_xdof_order, node2_remainder_xdof_order];
                

                //List<int> crosspoint_lagranges_X = new List<int>();
                int cross_point_remainder_xdof_order = subdRemainderDofs[model.GetNode(crossPointIds[0]), StructuralDof.TranslationX]; // TODO 1
                double krc_value2 = Krc[cross_point_remainder_xdof_order, cornerNode_xdof_order];
                //double krr_value2 = Krr[cross_point_remainder_xdof_order, node2_remainder_xdof_order];
                double krr_value3 = Krr[cross_point_remainder_xdof_order, cross_point_remainder_xdof_order];
                */
                #endregion
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                // Calculate the preconditioner before factorizing each subdomain's Krr.
                // The inter-subdomain stiffness distribution may have changed even if a subdomain's stiffness is the same.
                Logger.StartMeasuringTime();
                if (preconditioning.ReorderInternalDofsForFactorization) dofSeparator.ReorderInternalDofs(matrixManager);
                preconditioner = precondFactory.CreatePreconditioner(preconditioning, model, dofSeparator, lagrangesEnumerator,
                    matrixManager, stiffnessDistribution);
                Logger.LogCurrentTaskDuration("Calculating preconditioner");

                Logger.StartMeasuringTime();

                // Factorize each subdomain's Krr
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        //TODO: If I can reuse Krr, I can also reuse its factorization. Therefore this must be inPlace. In contrast, FETI-1 needs Kff intact for Stiffness distribution, in the current design).
                        Debug.WriteLine(msgHeader
                            + $"Inverting the remainder-remainder stiffness matrix of subdomain {subdomain.ID} in place.");
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).InvertKrr(true);
                        #region debug1
                        // TODOGer: edw ginetai extract Kcc, Krc kai Krr. en telei sto FetiDPsubdomainMatrixmanagerSkyline %115
                        #endregion
                    }
                }

                // Calculate FETI-DP coarse problem matrix
                matrixManager.CalcInverseCoarseProblemMatrix(CornerNodes);
                flexibility = new FetiDPFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, matrixManager);
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                isStiffnessModified = false;
            }

            // Calculate FETI-DP coarse problem rhs 
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

            // Solve interface problem
            Vector lagranges = interfaceProblemSolver.SolveInterfaceProblem(matrixManager, lagrangesEnumerator,
                flexibility, preconditioner, globalForcesNorm, Logger);
            Logger.LogCurrentTaskDuration("Solving interface problem");

            // Calculate the displacements of each subdomain
            Logger.StartMeasuringTime();
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);
            Logger.LogCurrentTaskDuration("Calculate displacements from lagrange multipliers");

            Logger.IncrementAnalysisStep();
        }

        public class Builder
        {
            private readonly IFetiDPMatrixManagerFactory matrixManagerFactory;

            public Builder(IFetiDPMatrixManagerFactory matrixManagerFactory)
            {
                this.matrixManagerFactory = matrixManagerFactory;
            }
            public ICrosspointStrategy CrosspointStrategy { get; set; } = new FullyRedundantConstraints();

            public IDofOrderer DofOrderer { get; set; } =
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public PcgSettings PcgSettings { get; set; } = new PcgSettings();

            public IFetiPreconditioningOperations Preconditioning { get; set; } = new DirichletPreconditioning();

            public StiffnessDistributionType StiffnessDistribution { get; set; } = StiffnessDistributionType.Homogeneous;

            public FetiDPSolverSerial Build(IModel model, ICornerNodeSelection cornerNodeSelection)
            {
                return new FetiDPSolverSerial(model, cornerNodeSelection, matrixManagerFactory, Preconditioning,
                    CrosspointStrategy, PcgSettings, StiffnessDistribution);
            }
        }
    }
}