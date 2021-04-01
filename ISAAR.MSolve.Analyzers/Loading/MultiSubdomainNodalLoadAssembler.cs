using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;

//TODO: I do not like the name "Assembler". I prefer "Applier"
//TODO: Right now each LoadApplier adds its contribution to a single vector (the rhs). Would it be better if it returned a vector 
//      and the analyzer added that vector to the rhs? That way contributions from some loads do not need to be recomputed 
//      everytime the rest do.
namespace ISAAR.MSolve.Analyzers.Loading
{
    public class MultiSubdomainNodalLoadAssembler : INodalLoadAssembler
    {
        private readonly INodalLoadDistributor distributor;

        public MultiSubdomainNodalLoadAssembler(INodalLoadDistributor distributor)
        {
            this.distributor = distributor;
        }

        public void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhsVector)
        {
            foreach (INodalLoad load in subdomain.EnumerateNodalLoads())
            {
                int idx = subdomain.FreeDofOrdering.FreeDofs[load.Node, load.DOF];
                double amount = distributor.ScaleNodalLoad(subdomain, load);
                rhsVector.Set(idx, amount);
            }
        }

        //public void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhs)
        //{
        //    //TODO: Scaling all dofs at once using SparseVector seems like a good idea, but SparseVector will be very slow for
        //    //      that operation.
        //    var loadVectorBuilder = new SortedDictionary<int, double>();
        //    foreach (INodalLoad load in subdomain.EnumerateNodalLoads())
        //    {
        //        INode node = load.Node;
        //        IDofType dof = load.DOF;

        //        // Add it to the vector builder
        //        int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dof];
        //        loadVectorBuilder[subdomainDofIdx] = load.Amount;
        //    }
        //    var loads = SparseVector.CreateFromDictionary(subdomain.FreeDofOrdering.NumFreeDofs, loadVectorBuilder); ;
        //    solver.ScaleNodalLoads(subdomain, loads);
        //    rhs.AddIntoThis(loads);
        //}
    }
}
