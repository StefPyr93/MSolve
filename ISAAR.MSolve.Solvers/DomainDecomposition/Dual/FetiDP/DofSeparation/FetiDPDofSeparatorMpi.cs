using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using MPI;

//TODO: Perhaps I should also find and expose the indices of boundary remainder and internal remainder dofs into the sequence 
//      of all free dofs of each subdomain
//TODO: Decide which of these data structures will be stored and which will be used ONCE to create all required mapping matrices.
//TODO: Perhaps the corner dof logic should be moved to another class.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPDofSeparatorMpi : IFetiDPDofSeparator
    {
        private readonly ICornerNodeSelection cornerNodeSelection;
        private readonly FetiDPGlobalDofSeparator globalDofs;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly ProcessDistribution procs;
        private readonly Dictionary<ISubdomain, FetiDPSubdomainDofSeparator> subdomainDofs;
        private readonly Dictionary<ISubdomain, DofTable> subdomainCornerDofOrderings_master = new Dictionary<ISubdomain, DofTable>();

        // These are defined per subdomain and are needed both in the corresponding process and in master.
        private Dictionary<ISubdomain, UnsignedBooleanMatrix> subdomainCornerBooleanMatrices_master;

        public FetiDPDofSeparatorMpi(ProcessDistribution processDistribution, IModel model,
            ICornerNodeSelection cornerNodeSelection)
        {
            this.procs = processDistribution;
            this.model = model;
            this.cornerNodeSelection = cornerNodeSelection;

            subdomainDofs = new Dictionary<ISubdomain, FetiDPSubdomainDofSeparator>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                subdomainDofs[subdomain] = new FetiDPSubdomainDofSeparator(subdomain);
            }
            if (procs.IsMasterProcess) globalDofs = new FetiDPGlobalDofSeparator(model);

            this.msgHeader = $"Process {processDistribution.OwnRank}, {this.GetType().Name}: ";
        }

        public Dictionary<INode, IDofType[]> GlobalBoundaryDofs
        {
            get
            {
                procs.CheckProcessIsMaster();
                return globalDofs.GlobalBoundaryDofs;
            }
        }

        public DofTable GlobalCornerDofOrdering
        {
            get
            {
                procs.CheckProcessIsMaster();
                return globalDofs.GlobalCornerDofOrdering;
            }
        }

        public int[] GlobalCornerToFreeDofMap
        {
            get
            {
                procs.CheckProcessIsMaster();
                return globalDofs.GlobalCornerToFreeDofMap;
            }
        }

        public int NumGlobalCornerDofs { get; private set; }

        public int[] GetBoundaryDofIndices(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].BoundaryDofIndices;
        }

        public (INode node, IDofType dofType)[] GetBoundaryDofs(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].BoundaryDofs;
        }

        public UnsignedBooleanMatrix GetCornerBooleanMatrix(ISubdomain subdomain)
        {
            if (procs.IsMasterProcess) return subdomainCornerBooleanMatrices_master[subdomain];
            else
            {
                procs.CheckProcessMatchesSubdomain(subdomain.ID);
                return subdomainDofs[subdomain].CornerBooleanMatrix;
            }
        }

        public int[] GetCornerDofIndices(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].CornerDofIndices;
        }

        public DofTable GetCornerDofOrdering(ISubdomain subdomain)
        {
            if (procs.IsMasterProcess) return subdomainCornerDofOrderings_master[subdomain];
            else
            {
                procs.CheckProcessMatchesSubdomain(subdomain.ID);
                return subdomainDofs[subdomain].CornerDofOrdering;
            }
        }

        public IReadOnlyList<(INode node, IDofType dofType)> GetCornerDofs(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].GetCornerDofs(cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain));
        }

        public DofTable GetRemainderDofOrdering(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].RemainderDofOrdering;
        }

        public int[] GetRemainderDofIndices(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].RemainderDofIndices;
        }

        public int[] GetInternalDofIndices(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return subdomainDofs[subdomain].InternalDofIndices;
        }

        public void ReorderInternalDofs(IFetiDPSeparatedDofReordering reordering)
        {
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                if (subdomain.ConnectivityModified)
                {
                    Debug.WriteLine(msgHeader + $"Reordering internal dofs of subdomain {subdomain.ID}.");
                    subdomainDofs[subdomain].ReorderInternalDofs(reordering.ReorderSubdomainInternalDofs(subdomain));
                }
            }
        }

        //TODO: This is too detailed for a coordinator class. All calls to global separator should be in 1 method. Ditto for subdomain separator. 
        public void SeparateDofs(IFetiDPSeparatedDofReordering reordering)
        {
            // Global dofs
            if (procs.IsMasterProcess)
            {
                globalDofs.DefineGlobalBoundaryDofs(cornerNodeSelection.GlobalCornerNodes);
                globalDofs.DefineGlobalCornerDofs(cornerNodeSelection.GlobalCornerNodes);
                NumGlobalCornerDofs = globalDofs.NumGlobalCornerDofs;
            }

            // Subdomain dofs
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                if (subdomain.ConnectivityModified)
                {
                    HashSet<INode> cornerNodes = cornerNodeSelection.GetCornerNodesOfSubdomain(subdomain);

                    Debug.WriteLine(msgHeader + $"Separating and ordering corner-remainder dofs of subdomain {s}");
                    subdomainDofs[subdomain].SeparateCornerRemainderDofs(cornerNodes);

                    Debug.WriteLine(msgHeader + $"Reordering internal dofs of subdomain {s}.");
                    subdomainDofs[subdomain].ReorderRemainderDofs(reordering.ReorderSubdomainRemainderDofs(subdomain));

                    Debug.WriteLine(msgHeader + $"Separating and ordering boundary-internal dofs of subdomain {s}");
                    subdomainDofs[subdomain].SeparateBoundaryInternalDofs(cornerNodes);
                }
            }

            // Reorder global corner dofs
            GatherCornerDofOrderingsFromSubdomains();
            if (procs.IsMasterProcess) globalDofs.ReorderGlobalCornerDofs(reordering.ReorderGlobalCornerDofs());

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
            if (procs.IsMasterProcess)
            {
                subdomainCornerBooleanMatrices_master = globalDofs.CalcCornerMappingMatrices(subdomainCornerDofOrderings_master);
            }
            ScatterCornerBooleanMatricesToSubdomains();
        }

        //TODO: The solver defines which subdomains are modified, but how? 
        //      Perhaps the object that provides corner nodes should also inform if they have changed.
        //TODO: There are alternative ways to get these in master process. E.g. The master process could redo the work required 
        //      to create it. Not sure if it will be slower. Should I make this a strategy? If you think about it, this is 
        //      probably the only thing that is fundamentally different between serial and MPI implementations. 
        //TODO: Add optimization in case all subdomains are modified (e.g. at the start of an analysis). In this case, use
        //      MPI gather, since it is faster than send/receive.
        private void GatherCornerDofOrderingsFromSubdomains()
        {
            var tableSerializer = new DofTableSerializer(model.DofSerializer);
            var activeSubdomains = new ActiveSubdomains(procs, s => model.GetSubdomain(s).ConnectivityModified);
            
            // Prepare data in each process
            var processOrderings = new Dictionary<int, DofTable>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                processOrderings[s] = subdomainDofs[subdomain].CornerDofOrdering;
            }

            // Gather data in master
            var transferrer = new TransferrerPerSubdomain(procs);
            GetArrayLengthOfPackedData<DofTable> getPackedDataLength = (s, table) => tableSerializer.CalcPackedLength(table);
            PackSubdomainDataIntoArray<DofTable, int> packData =
                (s, table, buffer, offset) => tableSerializer.PackTableIntoArray(table, buffer, offset);
            UnpackSubdomainDataFromArray<DofTable, int> unpackData =
                (s, buffer, start, end) => tableSerializer.UnpackTableFromArray(buffer, start, end, model.GetNode);
            Dictionary<int, DofTable> allOrderings = transferrer.GatherFromSomeSubdomainsPacked(processOrderings, 
                getPackedDataLength, packData, unpackData, activeSubdomains);

            // Store the received data in master
            //if (procs.IsMasterProcess) subdomainCornerDofOrderings_master = allOrderings.ChangeKey(model);
            if (procs.IsMasterProcess)
            {
                foreach (var idOrderingPair in allOrderings)
                {
                    ISubdomain subdomain = model.GetSubdomain(idOrderingPair.Key);
                    subdomainCornerDofOrderings_master[subdomain] = idOrderingPair.Value;
                }
            }
        }

        /// <summary>
        /// This should not be called, unless at least one corner dof of the whole model is modified.
        /// </summary>
        /// <param name="subdomain"></param>
        private void ScatterCornerBooleanMatricesToSubdomains()
        {
            // Scatter the matrices from master process
            var transferrer = new TransferrerPerSubdomain(procs);
            Dictionary<int, UnsignedBooleanMatrix> allMatricesBc = null;
            if (procs.IsMasterProcess) //TODO: Perhaps I should make the Dictionary have int as keys.
            {
                allMatricesBc = subdomainCornerBooleanMatrices_master.ChangeKey();
            }
            Dictionary<int, UnsignedBooleanMatrix> processMatricesBc = transferrer.ScatterToAllSubdomains(allMatricesBc);

            // Store them in other processes
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                var temp = subdomainDofs[subdomain];
                var Bc = processMatricesBc[s];
                temp.SetCornerBooleanMatrix(Bc, this);
                NumGlobalCornerDofs = Bc.NumColumns;
            }
        }
    }
}
