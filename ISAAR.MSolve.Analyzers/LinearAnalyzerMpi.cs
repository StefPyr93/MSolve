using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzerMpi : IChildAnalyzer
    {
        private readonly IModel model;
        private readonly ProcessDistribution procs;
        private readonly IAnalyzerProvider provider;
        private readonly ISolverMpi solver;

        public LinearAnalyzerMpi(ProcessDistribution processDistribution, IModel model, ISolverMpi solver, 
            IAnalyzerProvider provider)
        {
            this.procs = processDistribution;
            this.model = model;
            this.solver = solver;
            this.provider = provider;
        }

        public ILogFactory LogFactory { get; set; }
        public IAnalyzerLog[] Logs { get; set; }

        public IParentAnalyzer ParentAnalyzer { get; set; }

        Dictionary<int, IAnalyzerLog[]> IAnalyzer.Logs
        {
            get
            {
                //TODO: I should probably gather logs from processes here.
                throw new NotImplementedException();
            }
        }

        public void BuildMatrices()
        {
            if (ParentAnalyzer == null) throw new InvalidOperationException("This linear analyzer has no parent.");

            ParentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            //InitializeLogs();
            //solver.Initialize(); //TODO: Using this needs refactoring
        }

        public void Solve()
        {
            DateTime start = DateTime.Now;
            AddEquivalentNodalLoadsToRHS(); //TODO: The initial rhs (from other loads) should also be built by the analyzer instead of the model.
            solver.Solve();
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        private void AddEquivalentNodalLoadsToRHS()
        {
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                var rhs = solver.GetLinearSystem(subdomain).RhsVector;
                provider.DirichletLoadsAssembler.ApplyEquivalentNodalLoads(subdomain, rhs);
            }
        }

        private void InitializeLogs()
        {
            Logs = LogFactory.CreateLogs();
        }

        private void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                foreach (IAnalyzerLog log in Logs) log.StoreResults(start, end, solver.GetLinearSystem(subdomain).Solution);
            }
        }
    }
}
