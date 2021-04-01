using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.Loading
{
    public class SingleSubdomainNodalLoadAssembler : INodalLoadAssembler
    {
        public void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhsVector)
        {
            foreach (INodalLoad load in subdomain.EnumerateNodalLoads())
            {
                int idx = subdomain.FreeDofOrdering.FreeDofs[load.Node, load.DOF];
                rhsVector.Set(idx, load.Amount);
            }
        }
    }
}
