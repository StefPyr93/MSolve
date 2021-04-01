using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Right now each it adds its contribution to a single vector (the rhs). Would it be better if it returned a vector 
//      and the analyzer added that vector to the rhs? That way contributions from some loads do not need to be recomputed 
//      everytime the rest do. Then I would need an efficient and safe ZeroVector : IVectorView, rather than doing hacks
//      with SparseVector.
//TODO: I do not like the name "Assembler". "Applier" perhaps
namespace ISAAR.MSolve.Analyzers.Loading
{
    public interface INodalLoadAssembler
    {
        void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhsVector);
    }
}
