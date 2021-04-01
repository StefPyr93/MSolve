using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Implementations of this should probably have state (e.g. the corner nodes) which can be updated or not. The solver should 
//      query them if the corner nodes have changed, so that it will know that redefining dofs and lagrange multipliers must be 
//      performed. Accessing the corner nodes should be O(1)
//TODO: Clients should not have to implement both serial and MPI versions of their ICornerNodeSelection implementations.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public interface ICornerNodeSelection
    {
        HashSet<INode> GlobalCornerNodes { get; }

        HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain);

        void Update(); // TODO: This should notify FETI-DP solver when then corner nodes are actually updated.
    }
}
