using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    public class LinearAnalyzer : IChildAnalyzer
    {
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly IModel model;
        private readonly IAnalyzerProvider provider;
        private readonly ISolver solver;

        public LinearAnalyzer(IModel model, ISolver solver, IAnalyzerProvider provider)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
        }

        public Dictionary<int, ILogFactory> LogFactories { get; } = new Dictionary<int, ILogFactory>();
        public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

        public IParentAnalyzer ParentAnalyzer { get; set; }

        public void BuildMatrices()
        {
            if (ParentAnalyzer == null) throw new InvalidOperationException("This linear analyzer has no parent.");

            ParentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        public void Initialize(bool isFirstAnalysis)
        {
            InitializeLogs();
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
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                provider.DirichletLoadsAssembler.ApplyEquivalentNodalLoads(linearSystem.Subdomain, linearSystem.RhsVector);
            }
        }

        private void InitializeLogs()
        {
            Logs.Clear();
            foreach (int id in LogFactories.Keys) Logs.Add(id, LogFactories[id].CreateLogs());
        }

        private void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in Logs.Keys)
                foreach (var l in Logs[id])
                    l.StoreResults(start, end, linearSystems[id].Solution);
        }
    }
}
