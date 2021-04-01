using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;

//TODO: Remove code duplication between this and Feti1DofSeparator
//TODO: Perhaps I should also find and expose the indices of boundary remainder and internal remainder dofs into the sequence 
//      of all free dofs of each subdomain
//TODO: Decide which of these data structures will be stored and which will be used ONCE to create all required mapping matrices.
//TODO: Perhaps the corner dof logic should be moved to another class.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPDofSeparatorSerial : IFetiDPDofSeparator
    {
        private readonly ICornerNodeSelection cornerNodeSelection;
        private readonly FetiDPGlobalDofSeparator globalDofs;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly Dictionary<ISubdomain, FetiDPSubdomainDofSeparator> subdomainDofs;

        private Dictionary<ISubdomain, DofTable> subdomainCornerDofOrderings;

        public FetiDPDofSeparatorSerial(IModel model, ICornerNodeSelection cornerNodeSelection)
        {
            this.model = model;
            this.cornerNodeSelection = cornerNodeSelection;

            this.globalDofs = new FetiDPGlobalDofSeparator(model);
            this.subdomainDofs = new Dictionary<ISubdomain, FetiDPSubdomainDofSeparator>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                this.subdomainDofs[subdomain] = new FetiDPSubdomainDofSeparator(subdomain);
            }

            this.msgHeader = $"{this.GetType().Name}: ";
        }

        public Dictionary<INode, IDofType[]> GlobalBoundaryDofs => globalDofs.GlobalBoundaryDofs;
        public DofTable GlobalCornerDofOrdering => globalDofs.GlobalCornerDofOrdering;
        public int[] GlobalCornerToFreeDofMap => globalDofs.GlobalCornerToFreeDofMap;
        public int NumGlobalCornerDofs => globalDofs.NumGlobalCornerDofs;

        public int[] GetBoundaryDofIndices(ISubdomain subdomain) => subdomainDofs[subdomain].BoundaryDofIndices;
        public (INode node, IDofType dofType)[] GetBoundaryDofs(ISubdomain subdomain) => subdomainDofs[subdomain].BoundaryDofs;
        public UnsignedBooleanMatrix GetCornerBooleanMatrix(ISubdomain subdomain) => subdomainDofs[subdomain].CornerBooleanMatrix;
        public int[] GetCornerDofIndices(ISubdomain subdomain) => subdomainDofs[subdomain].CornerDofIndices;
        public DofTable GetCornerDofOrdering(ISubdomain subdomain) => subdomainDofs[subdomain].CornerDofOrdering;

        public IReadOnlyList<(INode node, IDofType dofType)> GetCornerDofs(ISubdomain subdomain)
            => subdomainDofs[subdomain].GetCornerDofs(cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain));

        public int[] GetInternalDofIndices(ISubdomain subdomain) => subdomainDofs[subdomain].InternalDofIndices;
        public DofTable GetRemainderDofOrdering(ISubdomain subdomain) => subdomainDofs[subdomain].RemainderDofOrdering;
        public int[] GetRemainderDofIndices(ISubdomain subdomain) => subdomainDofs[subdomain].RemainderDofIndices;

        public void ReorderInternalDofs(IFetiDPSeparatedDofReordering reordering)
        {
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.ConnectivityModified)
                {
                    Debug.WriteLine(msgHeader + $"Reordering internal dofs of subdomain {subdomain.ID}.");
                    subdomainDofs[subdomain].ReorderInternalDofs(reordering.ReorderSubdomainInternalDofs(subdomain));
                }
            }
        }

        public void SeparateDofs(IFetiDPSeparatedDofReordering reordering)
        {
            // Global dofs
            globalDofs.DefineGlobalBoundaryDofs(cornerNodeSelection.GlobalCornerNodes);
            globalDofs.DefineGlobalCornerDofs(cornerNodeSelection.GlobalCornerNodes);

            // Subdomain dofs
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                HashSet<INode> cornerNodes = cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain);
                if (subdomain.ConnectivityModified)
                {
                    int s = subdomain.ID;
                    Debug.WriteLine(msgHeader + $"Separating and ordering corner-remainder dofs of subdomain {s}");
                    subdomainDofs[subdomain].SeparateCornerRemainderDofs(cornerNodes);

                    Debug.WriteLine(msgHeader + $"Reordering internal dofs of subdomain {s}.");
                    subdomainDofs[subdomain].ReorderRemainderDofs(reordering.ReorderSubdomainRemainderDofs(subdomain));

                    Debug.WriteLine(msgHeader + $"Separating and ordering boundary-internal dofs of subdomain {s}");
                    subdomainDofs[subdomain].SeparateBoundaryInternalDofs(cornerNodes);
                }
            }

            // Reorder global corner dofs
            GatherSubdomainCornerDofOrderings();
            globalDofs.ReorderGlobalCornerDofs(reordering.ReorderGlobalCornerDofs());

            // Subdomain - global mappings
            CalcCornerMappingMatrices();
        }

        /// <summary>
        /// Bc unsigned boolean matrices that map global to subdomain corner dofs. This method must be called after 
        /// <see cref="DefineGlobalCornerDofs(Dictionary{int, HashSet{INode}})"/> and after the reordering of the global corner 
        /// dofs has been calculated.
        /// </summary>
        private void CalcCornerMappingMatrices()
        {
            // Create the corner mapping matrices
            Dictionary<ISubdomain, UnsignedBooleanMatrix> subdomainCornerBooleanMatrices
                = globalDofs.CalcCornerMappingMatrices(subdomainCornerDofOrderings);

            // Assign them to the subdomain dof separators as well
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainDofs[subdomain].SetCornerBooleanMatrix(subdomainCornerBooleanMatrices[subdomain], this);
            }
        }

        private void GatherSubdomainCornerDofOrderings()
        {
            subdomainCornerDofOrderings = new Dictionary<ISubdomain, DofTable>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainCornerDofOrderings[subdomain] = subdomainDofs[subdomain].CornerDofOrdering;
            }
        }
    }
}
