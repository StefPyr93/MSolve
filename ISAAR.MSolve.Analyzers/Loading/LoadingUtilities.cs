using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Analyzers.Loading
{
    public static class LoadingUtilities
    {
        //TODO: This method adds contributions into ISubdomain.Forces and is called by parent analyzers (Static, Newmark, 
        //      ThermalDynamic). Conversely DirichletEquivalentLoadsAssembler is applied to ILinearSystem.RhsVector and called
        //      by child analyzers (Linear, LoadControl, DisplacementControl). They should be used similarly.
        public static void ApplyNodalLoads(IModel model, ISolver solver)
        {
            INodalLoadAssembler loadAssembler;
            if (model.NumSubdomains == 1) loadAssembler = new SingleSubdomainNodalLoadAssembler();
            else loadAssembler = new MultiSubdomainNodalLoadAssembler(solver.NodalLoadDistributor);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                loadAssembler.ApplyEquivalentNodalLoads(subdomain, subdomain.Forces);
            }
        }

        //TODO: This method adds contributions into ISubdomain.Forces and is called by parent analyzers (Static, Newmark, 
        //      ThermalDynamic). Conversely DirichletEquivalentLoadsAssembler is applied to ILinearSystem.RhsVector and called
        //      by child analyzers (Linear, LoadControl, DisplacementControl). They should be used similarly.
        public static void ApplyNodalLoads(IModel model, ISolverMpi solver)
        {
            var loadAssembler = new MultiSubdomainNodalLoadAssembler(solver.NodalLoadDistributor);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                loadAssembler.ApplyEquivalentNodalLoads(subdomain, subdomain.Forces);
            }
        }

        public static void ApplyNodalLoadsMpi(ProcessDistribution processDistribution, IModel model, ISolverMpi solver)
        {
            var loadAssembler = new MultiSubdomainNodalLoadAssembler(solver.NodalLoadDistributor);
            foreach (int s in processDistribution.GetSubdomainIDsOfProcess(processDistribution.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                loadAssembler.ApplyEquivalentNodalLoads(subdomain, subdomain.Forces);
            }
        }
    }
}
