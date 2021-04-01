using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Perhaps it should be joined with INodalLoadAssembler.
namespace ISAAR.MSolve.Analyzers.Loading
{
    public interface IDirichletEquivalentLoadsAssembler
    {
        void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhs);
    }
}
