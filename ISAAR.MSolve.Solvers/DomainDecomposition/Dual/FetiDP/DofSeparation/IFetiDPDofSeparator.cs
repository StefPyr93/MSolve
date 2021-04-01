using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;

//TODO: Perhaps there should only be 1 method for building the items that are 100% necessary, instead of the solver having to call
//      many other ones in the correct order. Only the reordering cannot be performed by this dof separators, so it should be 
//      injected as a delegate or as an interface (e.g. IMatrixManager or something above for looser coupling).
//TODO: Perhaps corner nodes should be accessed by ICornerNodesSelection. Otherwise MPI processes other then master would 
//      have to pass null for global ones and this would be done by yet another if (IsMasterProcess) ...
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public interface IFetiDPDofSeparator : IDofSeparator
    {
        /// <summary>
        /// Dof ordering for corner dofs of the model: Each (INode, IDofType) pair is associated with the index of that dof into 
        /// a vector corresponding to all corner dofs of the model.
        /// </summary>
        DofTable GlobalCornerDofOrdering { get; }

        /// <summary>
        /// If Xf is a vector with all free dofs of the model and Xc is a vector with all corner dofs of the model, then
        /// Xf[GlobalCornerToFreeDofMap[i]] = Xc[i].
        /// </summary>
        int[] GlobalCornerToFreeDofMap { get; }

        /// <summary>
        /// The number of corner dofs of the model.
        /// </summary>
        int NumGlobalCornerDofs { get; }

        /// <summary>
        /// Corner dofs of each subdomain. They have the same order as <see cref="GetCornerDofIndices(ISubdomain)"/>.
        /// </summary>
        IReadOnlyList<(INode node, IDofType dofType)> GetCornerDofs(ISubdomain subdomain);

        UnsignedBooleanMatrix GetCornerBooleanMatrix(ISubdomain subdomain);

        /// <summary>
        /// Indices of (boundary) corner dofs into the sequence of all free dofs of a subdomain.
        /// </summary>
        int[] GetCornerDofIndices(ISubdomain subdomain);

        /// <summary>
        /// Dof ordering for corner dofs of a subdomain: Each (INode, IDofType) pair of the subdomain is associated with the   
        /// index of that dof into a vector corresponding to corner dofs of that subdomain.
        /// </summary>
        DofTable GetCornerDofOrdering(ISubdomain subdomain);

        /// <summary>
        /// Dof ordering for remainder (boundary and internal) dofs of a subdomain: Each (INode, IDofType) pair of the 
        /// subdomain is associated with the index of that dof into a vector corresponding to remainder dofs of that subdomain.
        /// </summary>
        DofTable GetRemainderDofOrdering(ISubdomain subdomain);

        /// <summary>
        /// Indices of remainder (boundary and internal) dofs into the sequence of all free dofs of a subdomain.
        /// </summary>
        int[] GetRemainderDofIndices(ISubdomain subdomain);

        void SeparateDofs(IFetiDPSeparatedDofReordering reordering);

        void ReorderInternalDofs(IFetiDPSeparatedDofReordering reordering);
    }
}
