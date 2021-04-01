using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using MPI;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;

// TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip can be found
namespace ISAAR.MSolve.XFEM.Analyzers
{
    /// <summary>
    /// Implements crack propagation under static loading with linear material behavior. Based on Linear Elastic Fracture 
    /// Mechanics. Appropriate for brittle materials or fatigue crack propagation analysis. For now, it only works with XFEM.
    /// </summary>
    public class QuasiStaticCrackPropagationAnalyzerMpiRedundnat //: IAnalyzer
    {
        private readonly ICrackDescriptionMpi crack;
        private readonly double fractureToughness;
        private readonly int maxIterations;
        private readonly IXModelMpi model;
        private readonly TipAdaptivePartitioner partitioner; //TODO: Refactor its injection and usage
        private readonly ProcessDistribution procs;
        private readonly bool reanalysis;

        //private readonly IStaticProvider problem; //TODO: refactor and use this instead
        private readonly ElementStructuralStiffnessProvider problem = new ElementStructuralStiffnessProvider();
        private readonly DirichletEquivalentLoadsAssembler loadsAssembler; 

        private readonly ISolverMpi solver;
        private HashSet<ISubdomain> newTipEnrichedSubdomains;
        private CrackPropagationTermination termination;

        public QuasiStaticCrackPropagationAnalyzerMpiRedundnat(ProcessDistribution processDistribution, IXModelMpi model, 
            ISolverMpi solver, /*IStaticProvider problem,*/ ICrackDescriptionMpi crack, double fractureToughness, 
            int maxIterations, bool reanalysis = true, TipAdaptivePartitioner partitioner = null)
        {
            this.procs = processDistribution;
            this.model = model;
            this.solver = solver;
            //this.problem = problem;
            this.crack = crack;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
            this.reanalysis = reanalysis;
            this.partitioner = partitioner;

            //TODO: Refactor problem structural and remove the next
            problem = new ElementStructuralStiffnessProvider();
            loadsAssembler = new DirichletEquivalentLoadsAssembler(problem); ;
        }

        public IDomainDecompositionLogger DDLogger { get; set; }
        public CrackPropagationTermination Termination => termination;

        public void Initialize(bool isFirstAnalysis = true)
        {
            // The order in which the next initializations happen is very important.
            if (isFirstAnalysis) model.ConnectDataStructures();
            //model.ScatterSubdomains(); //This is not needed, since each process stores the whole model

            //solver.Initialize(); //TODO: not sure about this one.
        }

        /// <summary>
        /// Returns the crack path after repeatedly executing: XFEM analysis, SIF calculation, crack propagation
        /// </summary>
        /// <returns></returns>
        public void Analyze()
        {
            int analysisStep;
            for (analysisStep = 0; analysisStep < maxIterations; ++analysisStep)
            {
                procs.Communicator.Barrier();
                if (procs.IsMasterProcess)
                {
                    Debug.WriteLine($"Process {procs.OwnRank}: Crack propagation step {analysisStep}");
                    Console.WriteLine($"Process {procs.OwnRank}: Crack propagation step {analysisStep}. Tip ({crack.CrackTips[0].X}, {crack.CrackTips[0].Y})");
                }

                // Apply the enrichments due to the updated crack
                //Console.WriteLine($"Process {procs.OwnRank}: Update enrichment data");
                crack.UpdateEnrichments();

                // Update the subdomains due to the new enrichments and possible repartition the mesh
                //Console.WriteLine($"Process {procs.OwnRank}: Update mesh partition");
                UpdateSubdomains();

                //TODO: This is not needed since each process stores and updates the whole crack.
                //// Scatter the crack and enrichment data to other processes
                ////Console.WriteLine($"Process {procs.OwnRank}: Scattering crack data");
                //crack.ScatterCrackData(model); 

                // Order and count dofs
                //Console.WriteLine($"Process {procs.OwnRank}: Ordering dofs and crack data");
                solver.OrderDofs(false);
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    ISubdomain subdomain = model.GetSubdomain(s);
                    ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }

                // Create the stiffness matrix and then the forces vector
                //Console.WriteLine($"Process {procs.OwnRank}: Calculating matrix and rhs");
                //problem.ClearMatrices();
                solver.BuildGlobalMatrix(problem);
                //PrintKff(procs, linearSystem);
                model.ApplyLoads();
                LoadingUtilities.ApplyNodalLoadsMpi(procs, model, solver);
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    ISubdomain subdomain = model.GetSubdomain(s);
                    ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                    loadsAssembler.ApplyEquivalentNodalLoads(subdomain, linearSystem.RhsVector);
                }
                
                // Plot domain decomposition data, if necessary
                if (procs.IsMasterProcess)
                {
                    if (DDLogger != null) DDLogger.PlotSubdomains(model);
                }

                // Solve the linear system
                //Console.WriteLine($"Process {procs.OwnRank}: Solving linear system(s)");
                solver.Solve();

                //// Output field data
                //if (fieldOutput != null)
                //{
                //    fieldOutput.WriteOutputData(solver.DofOrderer, freeDisplacements, constrainedDisplacements, iteration);
                //}

                // Let the crack propagate
                //Console.WriteLine($"Process {procs.OwnRank}: Propagating the crack.");
                Dictionary<int, Vector> freeDisplacements = GatherDisplacementsToMaster();
                GatherSubdomainFreeDofOrderingsToMaster();
                crack.Propagate(freeDisplacements);

                // Check convergence 
                //Console.WriteLine($"Process {procs.OwnRank}: Checking convergence.");
                bool mustTerminate = MustTerminate(analysisStep);
                if (mustTerminate) throw new Exception("Early termination");
            }
            termination = CrackPropagationTermination.RequiredIterationsWereCompleted;
        }

        private static void PrintKff(ProcessDistribution procs, ILinearSystem linearSystem)
        {
            procs.Communicator.Barrier();
            string msg;
            try
            {
                int m = linearSystem.Matrix.NumRows;
                int n = linearSystem.Matrix.NumColumns;
                double norm = LinearAlgebra.Reduction.Reductions.Norm2(linearSystem.Matrix);
                msg = $"Subdomain {linearSystem.Subdomain.ID}: Kff ({m} x {n}), norm(Kff) = {norm}";
            }
            catch (Exception)
            {
                msg = $"Subdomain {linearSystem.Subdomain.ID}: Kff is the same as previous propagation step, but has been overwritten";
            }
            
            MpiUtilities.DoInTurn(procs.Communicator, () => Console.WriteLine(msg));
        }

        // TODO: Abstract this and add Tanaka_1974 approach
        private double CalculateEquivalentSIF(double sifMode1, double sifMode2)
        {
            return Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
        }

        private HashSet<ISubdomain> FindSubdomainsWithNewHeavisideEnrichedNodes()
        {
            var newHeavisideEnrichedSubdomains = new HashSet<ISubdomain>();
            foreach (ISet<XNode> heavisideNodes in crack.CrackBodyNodesNew.Values)
            {
                foreach (XNode node in heavisideNodes)
                {
                    newHeavisideEnrichedSubdomains.UnionWith(node.SubdomainsDictionary.Values);
                }
            }
            return newHeavisideEnrichedSubdomains;
        }

        private HashSet<ISubdomain> FindSubdomainsWithNewTipEnrichedNodes()
        {
            var newTipEnrichedSubdomains = new HashSet<ISubdomain>();
            foreach (ISet<XNode> tipNodes in crack.CrackTipNodesNew.Values)
            {
                foreach (XNode node in tipNodes)
                {
                    newTipEnrichedSubdomains.UnionWith(node.SubdomainsDictionary.Values);
                }
            }
            return newTipEnrichedSubdomains;
        }

        private Dictionary<int, Vector> GatherDisplacementsToMaster()
        {
            var transferrer = new TransferrerPerSubdomain(procs);
            var processUf = new Dictionary<int, Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                processUf[s] = (Vector)solver.GetLinearSystem(subdomain).Solution;
            }
            return transferrer.GatherFromAllSubdomains(processUf);
        }

        private void GatherSubdomainFreeDofOrderingsToMaster() //TODO: This should not be necessary
        {
            var globalDofOrdering = (GlobalFreeDofOrderingMpi)model.GlobalDofOrdering;
            //globalDofOrdering.GatherSubdomainDofOrderings(); // This should already have been called, during calculating the normalization of redisual in PCG.
            if (procs.IsMasterProcess)
            {
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    subdomain.FreeDofOrdering = model.GlobalDofOrdering.GetSubdomainDofOrdering(subdomain);
                }
            }
        }

        private bool MustTerminate(int analysisStep)
        {
            bool mustTerminate = false;
            if (procs.IsMasterProcess) mustTerminate = MustTerminate_master(analysisStep);
            procs.Communicator.Broadcast(ref mustTerminate, procs.MasterProcess);
            //procs.Communicator.Broadcast(ref termination, procs.MasterProcess); //TODO: This needs serialization of the enum and might not be necessary
            return mustTerminate;
        }

        private bool MustTerminate_master(int analysisStep)
        {
            // Check convergence 
            //TODO: Perhaps this should be done by the crack geometry or the Propagator itself and handled via exceptions 
            foreach (var tipPropagator in crack.CrackTipPropagators)
            {
                double sifEffective = CalculateEquivalentSIF(tipPropagator.Value.Logger.SIFsMode1[analysisStep],
                    tipPropagator.Value.Logger.SIFsMode2[analysisStep]);
                //Console.WriteLine("Keff = " + sifEffective);
                if (sifEffective >= fractureToughness)
                {
                    termination = CrackPropagationTermination.FractureToughnessIsExceeded;
                    return true;
                }
                if (!model.Boundary.IsInside(tipPropagator.Key))
                {
                    termination = CrackPropagationTermination.CrackExitsDomainBoundary;
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// Update the mesh partitioning and identify unmodified subdomains to avoid fully processing them again.
        /// </summary>
        /// <param name="repartitionedSubdomains"></param>
        private void UpdateSubdomains()
        {
            HashSet<ISubdomain> repartitionedSubdomains = null;

            if (model.NumSubdomains == 1) throw new InvalidOperationException("There must be >1 subdomains in a MPI environment");

            if (newTipEnrichedSubdomains == null) 
            {
                // First analysis step: All subdomains must be fully processed.
                if (partitioner != null) repartitionedSubdomains = partitioner.UpdateSubdomains();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    subdomain.ConnectivityModified = true;
                    subdomain.StiffnessModified = true;
                }

                // Prepare for the next analysis step
                newTipEnrichedSubdomains = FindSubdomainsWithNewTipEnrichedNodes();
            }
            else
            {
                // Update the mesh partitioning, if necessary
                HashSet<ISubdomain> modifiedSubdomains = new HashSet<ISubdomain>();
                if (partitioner != null)
                {
                    repartitionedSubdomains = partitioner.UpdateSubdomains();
                    modifiedSubdomains.UnionWith(repartitionedSubdomains);
                }

                if (reanalysis)
                {
                    // Set them all to unmodified.
                    foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                    {
                        subdomain.ConnectivityModified = false;
                        subdomain.StiffnessModified = false;
                    }

                    // The modified subdomains are the ones containing nodes enriched with tip or Heaviside functions during the  
                    // current analysis step. Also the ones that had tip enriched nodes in the previous step.
                    modifiedSubdomains.UnionWith(newTipEnrichedSubdomains);
                    newTipEnrichedSubdomains = FindSubdomainsWithNewTipEnrichedNodes(); // Prepare for the next analysis step
                    modifiedSubdomains.UnionWith(newTipEnrichedSubdomains);
                    HashSet<ISubdomain> newHeavisideSubdomains = FindSubdomainsWithNewHeavisideEnrichedNodes();
                    modifiedSubdomains.UnionWith(newHeavisideSubdomains);

                    foreach (ISubdomain subdomain in modifiedSubdomains)
                    {
                        subdomain.ConnectivityModified = true;
                        subdomain.StiffnessModified = true;
                    }
                }
            }
        }


        /// <summary>
        /// For debugging purposes
        /// </summary>
        /// <param name="subdomains"></param>
        /// <param name="name"></param>
        private void PrintSubdomainSubset(IEnumerable<ISubdomain> subdomains, string name)
        {
            var msg = new System.Text.StringBuilder($"Process {procs.OwnRank}: {name} subdomains = ");
            foreach (ISubdomain sub in subdomains) msg.Append(sub.ID + " ");
            Console.WriteLine(msg);
        }
    }
}
