using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Solvers;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Analyzers.ObjectManagers;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public class LoadControlAnalyzerDevelop3 : IChildAnalyzer //: NonLinearAnalyzerBaseDevelop2
    {

        IMaterialManager materialManager;

        public LoadControlAnalyzerDevelop3(IModel model, ISolver solver, INonLinearProvider provider,
            int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance, IMaterialManager materialManager)
        {
            this.model = model;
            this.solver = solver;
            this.provider = provider;
            this.subdomainUpdaters = CreateDefaultSubdomainUpdaters();
            this.linearSystems = solver.LinearSystems;
            this.numIncrements = numIncrements;
            this.maxIterationsPerIncrement = maxIterationsPerIncrement;
            this.numIterationsForMatrixRebuild = numIterationsForMatrixRebuild;
            this.residualTolerance = residualTolerance;
            this.materialManager = materialManager;
            //SubdomainUpdaters = CreateDefaultSubdomainUpdaters();
        }
        
        public void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();//TODOMaria this divides the externally applied load by the number of increments and scatters it to all subdomains and stores it in the class subdomain dictionary and total external load vector
            for (int increment = 0; increment < numIncrements; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();//TODOMaria this sets du to 0
                UpdateRhs(increment);//TODOMaria this copies the residuals stored in the class dictionary to the subdomains

                double firstError = 0;
                int iteration = 0;
                for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
                {
                    if (iteration == maxIterationsPerIncrement - 1) return;
                    if (Double.IsNaN(errorNorm)) return;
                    solver.Solve();
                    //double rhsNormIt = solver.LinearSystems.First().Value.RhsVector.Norm2();
                    //double xNormIt = solver.LinearSystems.First().Value.Solution.Norm2();


                    //Dictionary<int, IVector> internalRhsVectors = CalculateInternalRhs(increment, iteration);
                    CalculateInternalRhsStressesOnly(increment, iteration);
                    materialManager.UpdateMaterials();
                    Dictionary<int, IVector> internalRhsVectors = CalculateInternalRhsClaculateRhsOnly(increment, iteration);

                    double residualNormCurrent = UpdateResidualForcesAndNorm(increment, internalRhsVectors); // This also sets the rhs vectors in linear systems.
                    errorNorm = globalRhsNormInitial != 0 ? residualNormCurrent / globalRhsNormInitial : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)
                    //Console.WriteLine($"Increment {increment}, iteration {iteration}: norm2(error) = {errorNorm}");

                    if (iteration == 0) firstError = errorNorm;

                    if (TotalDisplacementsPerIterationLog != null) TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);

                    bool hasConverged;
                    hasConverged = errorNorm < residualTolerance;


                    if (hasConverged)
                    {
                        foreach (var subdomainLogPair in IncrementalLogs)
                        {
                            int subdomainID = subdomainLogPair.Key;
                            TotalLoadsDisplacementsPerIncrementLog log = subdomainLogPair.Value;
                            log.LogTotalDataForIncrement(increment, iteration, errorNorm,
                                uPlusdu[subdomainID], internalRhsVectors[subdomainID]);
                        }
                        break;
                    }

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((iteration + 1) % numIterationsForMatrixRebuild == 0) // Matrix rebuilding should be handled in another way. E.g. in modified NR, it must be done at each increment.
                    {
                        provider.Reset();
                        BuildMatrices();
                    }
                }
                //double rhsNormInc = solver.LinearSystems.First().Value.RhsVector.Norm2();
                //double xNormInc = solver.LinearSystems.First().Value.Solution.Norm2();
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);


                //SaveMaterialStateAndUpdateSolution();
                SaveSolution();
                materialManager.SaveState();


            }
            //            ClearMaterialStresses();

            // TODO: Logging should be done at each iteration. And it should be done using pull observers
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        #region fields from base copy
        protected readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        protected readonly int maxIterationsPerIncrement;
        protected readonly IModel model;
        protected readonly int numIncrements;
        protected readonly int numIterationsForMatrixRebuild;
        protected readonly INonLinearProvider provider;
        protected readonly double residualTolerance;
        protected readonly ISolver solver;
        protected readonly IReadOnlyDictionary<int, INonLinearSubdomainUpdaterDevelop> subdomainUpdaters;
        protected readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
        protected readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
        protected Vector globalRhs; //TODO: This was originally readonly 
        protected double globalRhsNormInitial; //TODO: This can probably be a local variable.
        protected INonLinearParentAnalyzer parentAnalyzer = null;

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get; } = new Dictionary<int, LinearAnalyzerLogFactory>();
        public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

        public TotalDisplacementsPerIterationLog TotalDisplacementsPerIterationLog { get; set; }
        public Dictionary<int, TotalLoadsDisplacementsPerIncrementLog> IncrementalLogs { get; }
            = new Dictionary<int, TotalLoadsDisplacementsPerIncrementLog>();

        public IParentAnalyzer ParentAnalyzer
        {
            get => parentAnalyzer;
            set => parentAnalyzer = (INonLinearParentAnalyzer)value; //TODO: remove this cast. Now it only serves as a check
        }
        #endregion

        #region methods from base copy

        public void BuildMatrices()
        {
            if (parentAnalyzer == null) throw new InvalidOperationException(
                "This Newton-Raphson nonlinear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            InitializeInternalVectors();
            //solver.Initialize(); //TODO: Using this needs refactoring
        }

        //TODO: Internal Rhs vectors are created and destroyed at each iteration. It would be more efficient to store them as
        //      vectors and then overwrite them.
        protected Dictionary<int, IVector> CalculateInternalRhs(int currentIncrement, int iteration)
        {
            var internalRhsVectors = new Dictionary<int, IVector>();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                if (currentIncrement == 0 && iteration == 0)
                {
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    du[id].Clear();
                    uPlusdu[id].Clear();
                    du[id].AddIntoThis(linearSystem.Solution);
                    uPlusdu[id].AddIntoThis(linearSystem.Solution);
                    du[id].SubtractIntoThis(u[id]);
                }
                else
                {

                    du[id].AddIntoThis(linearSystem.Solution);
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    uPlusdu[id].Clear();
                    uPlusdu[id].AddIntoThis(u[id]);
                    uPlusdu[id].AddIntoThis(du[id]);
                }
                //Vector<double> internalRhs = (Vector<double>)subdomain.GetRhsFromSolution(u[subdomain.ID], du[subdomain.ID]);

                //TODO: remove cast
                IVector internalRhs = subdomainUpdaters[id].GetRhsFromSolution(uPlusdu[id], du[id]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRhs(linearSystem.Subdomain, uPlusdu[id], internalRhs);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(linearSystem, uPlusdu[id]);
                    internalRhs.AddIntoThis(otherRhsComponents);//TODOMaria this does nothing for the static problem
                }

                internalRhsVectors.Add(id, internalRhs);
            }

            return internalRhsVectors;
        }

        protected void CalculateInternalRhsStressesOnly(int currentIncrement, int iteration)
        {
            var internalRhsVectors = new Dictionary<int, IVector>();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                if (currentIncrement == 0 && iteration == 0)
                {
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    du[id].Clear();
                    uPlusdu[id].Clear();
                    du[id].AddIntoThis(linearSystem.Solution);
                    uPlusdu[id].AddIntoThis(linearSystem.Solution);
                    du[id].SubtractIntoThis(u[id]);
                }
                else
                {

                    du[id].AddIntoThis(linearSystem.Solution);
                    //TODO: instead of Clear() and then AddIntoThis(), use only CopyFromVector()
                    uPlusdu[id].Clear();
                    uPlusdu[id].AddIntoThis(u[id]);
                    uPlusdu[id].AddIntoThis(du[id]);
                }
                //Vector<double> internalRhs = (Vector<double>)subdomain.GetRhsFromSolution(u[subdomain.ID], du[subdomain.ID]);

                subdomainUpdaters[id].CalculateStressesOnly(uPlusdu[id], du[id]);//TODOMaria this calculates the internal forces

            }

        }

        protected Dictionary<int, IVector> CalculateInternalRhsClaculateRhsOnly(int currentIncrement, int iteration)
        {
            var internalRhsVectors = new Dictionary<int, IVector>();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                //TODO: remove cast
                IVector internalRhs = subdomainUpdaters[id].CalculateRHSonly(uPlusdu[id], du[id]);//TODOMaria this calculates the internal forces
                provider.ProcessInternalRhs(linearSystem.Subdomain, uPlusdu[id], internalRhs);//TODOMaria this does nothing
                //(new Vector<double>(u[subdomain.ID] + du[subdomain.ID])).Data);

                if (parentAnalyzer != null)
                {
                    IVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(linearSystem, uPlusdu[id]);
                    internalRhs.AddIntoThis(otherRhsComponents);//TODOMaria this does nothing for the static problem
                }

                internalRhsVectors.Add(id, internalRhs);
            }

            return internalRhsVectors;
        }

        protected double UpdateResidualForcesAndNorm(int currentIncrement, Dictionary<int, IVector> internalRhs)
        {
            globalRhs.Clear();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                linearSystem.RhsVector.Clear(); //TODO: we can copy rhs[subdomain.ID] and then scale it instead of clearing and adding.

                // External forces = loadFactor * total external forces
                //TODO: the next line adds a vector to itself many times. This is called multiplication and is much faster.
                for (int j = 0; j <= currentIncrement; j++) linearSystem.RhsVector.AddIntoThis(rhs[id]);//TODOMaria this adds the external forces 

                // Residual forces = external - internal
                linearSystem.RhsVector.SubtractIntoThis(internalRhs[id]);

                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            return provider.CalculateRhsNorm(globalRhs);
        }

        protected void ClearIncrementalSolutionVector()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values) du[linearSystem.Subdomain.ID].Clear();
        }

        protected virtual void InitializeInternalVectors()//TODOMaria: this is probably where the initial internal nodal vector is calculated
        {
            globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            rhs.Clear();
            u.Clear();
            du.Clear();
            uPlusdu.Clear();

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)numIncrements);
                rhs.Add(id, r);
                u.Add(id, linearSystem.CreateZeroVector());
                du.Add(id, linearSystem.CreateZeroVector());
                uPlusdu.Add(id, linearSystem.CreateZeroVector());
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
        }

        protected void InitializeLogs()
        {
            Logs.Clear();
            foreach (int id in LogFactories.Keys) Logs.Add(id, LogFactories[id].CreateLogs());
            foreach (var log in IncrementalLogs.Values) log.Initialize();
        }

        protected void SaveMaterialStateAndUpdateSolution()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                subdomainUpdaters[id].UpdateState();
                u[id].AddIntoThis(du[id]);
            }
        }

        protected void SaveSolution()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                //subdomainUpdaters[id].UpdateState();
                u[id].AddIntoThis(du[id]);
            }
        }

        protected void SplitResidualForcesToSubdomains()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                linearSystem.RhsVector.Clear(); //TODO: why clear it if it is going to be overwritten immediately afterwards?
                model.GlobalDofOrdering.ExtractVectorSubdomainFromGlobal(linearSystem.Subdomain, globalRhs,
                    linearSystem.RhsVector);
            }
        }

        protected void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in Logs.Keys)
                foreach (var l in Logs[id])
                    l.StoreResults(start, end, u[id]);
        }

        protected void UpdateInternalVectors()//TODOMaria this is where I should add the calculation of the internal nodal force vector
        {
            globalRhs.Clear();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;

                //TODO: directly copy into linearSystem.RhsVector and then scale that.
                IVector r = linearSystem.RhsVector.Copy();
                r.ScaleIntoThis(1 / (double)numIncrements);
                rhs[id] = r;
                model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
            }
            globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
        }

        protected void UpdateRhs(int step)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);
                //linearSystem.RhsVector.Multiply(step + 1);
            }
        }

        private IReadOnlyDictionary<int, INonLinearSubdomainUpdaterDevelop> CreateDefaultSubdomainUpdaters()
        {
            int numSubdomains = model.NumSubdomains;
            var subdomainUpdaters = new Dictionary<int, INonLinearSubdomainUpdaterDevelop>(numSubdomains);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainUpdaters[subdomain.ID] = new NonLinearSubdomainUpdaterDevelop3(subdomain);

            }
            return subdomainUpdaters;

        }

        #endregion

        public class Builder
        {
            public Builder(IModel model, ISolver solver, INonLinearProvider provider, int numIncrements)
            {
                MaxIterationsPerIncrement = 1000;
                NumIterationsForMatrixRebuild = 1;
                ResidualTolerance = 1E-3;

                //TODO: this should belong to all (child) analyzer builders
                this.model = model;
                this.solver = solver;
                this.provider = provider;
                this.numIncrements = numIncrements;
                SubdomainUpdaters = CreateDefaultSubdomainUpdaters();

            }

            protected int maxIterationsPerIncrement = 1000;
            protected readonly IModel model;
            protected readonly int numIncrements;
            protected int numIterationsForMatrixRebuild = 1;
            protected readonly INonLinearProvider provider;
            protected double residualTolerance = 1e-8;
            protected readonly ISolver solver;

            //public LoadControlAnalyzerDevelop3 Build()
            //{
            //    return new LoadControlAnalyzerDevelop3(model, solver, provider, SubdomainUpdaters,
            //        numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
            //}

            public int MaxIterationsPerIncrement
            {
                get => maxIterationsPerIncrement;
                set
                {
                    if (value < 1) throw new ArgumentException("Max iterations per increment must be >= 1");
                    maxIterationsPerIncrement = value;
                }
            }

            public int NumIterationsForMatrixRebuild
            {
                get => numIterationsForMatrixRebuild;
                set
                {
                    if (value < 1) throw new ArgumentException("Iterations number for matrix rebuild must be >= 1");
                    numIterationsForMatrixRebuild = value;
                }
            }

            public double ResidualTolerance
            {
                get => residualTolerance;
                set
                {
                    if (value <= 0.0) throw new ArgumentException("Residual tolerance must be positive");
                    residualTolerance = value;
                }
            }

            public IReadOnlyDictionary<int, INonLinearSubdomainUpdaterDevelop> SubdomainUpdaters { get; set; }

            private IReadOnlyDictionary<int, INonLinearSubdomainUpdaterDevelop> CreateDefaultSubdomainUpdaters()
            {
                int numSubdomains = model.NumSubdomains;
                var subdomainUpdaters = new Dictionary<int, INonLinearSubdomainUpdaterDevelop>(numSubdomains);
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    subdomainUpdaters[subdomain.ID] = new NonLinearSubdomainUpdaterDevelop3(subdomain);

                }
                return subdomainUpdaters;

            }
        }
    }
}
