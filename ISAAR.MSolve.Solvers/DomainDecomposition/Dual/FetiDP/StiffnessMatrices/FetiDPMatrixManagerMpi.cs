using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Distributed.Transfer;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPMatrixManagerMpi : IFetiDPMatrixManager
    {
        private readonly IFetiDPGlobalMatrixManager matrixManagerGlobal_master;
        private readonly Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagerSubdomains;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly ProcessDistribution procs;

        public FetiDPMatrixManagerMpi(ProcessDistribution processDistribution, IModel model, IFetiDPDofSeparator dofSeparator,
            IFetiDPMatrixManagerFactory matrixManagerFactory)
        {
            this.procs = processDistribution;
            this.model = model;

            this.matrixManagerSubdomains = new Dictionary<int, IFetiDPSubdomainMatrixManager>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                this.matrixManagerSubdomains[s] = matrixManagerFactory.CreateSubdomainMatrixManager(subdomain, dofSeparator);
            }

            if (processDistribution.IsMasterProcess)
            {
                matrixManagerGlobal_master = matrixManagerFactory.CreateGlobalMatrixManager(model, dofSeparator);
            }
            this.msgHeader = $"Process {processDistribution.OwnRank}, {this.GetType().Name}: ";
        }

        public Vector CoarseProblemRhs
        {
            get
            {
                procs.CheckProcessIsMaster();
                return matrixManagerGlobal_master.CoarseProblemRhs;
            }
        }

        public void CalcCoarseProblemRhs()
        {
            // Calculate the subdomain FcStar in each process
            var processVectors = new Dictionary<int, Vector>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                // fcStar[s] = fbc[s] - Krc[s]^T * inv(Krr[s]) * fr[s] -> delegated to the SubdomainMatrixManager
                // globalFcStar = sum_over_s(Lc[s]^T * fcStar[s]) -> delegated to the GlobalMatrixManager
                matrixManagerSubdomains[s].CalcCoarseProblemRhsSubvectors();
                processVectors[s] = matrixManagerSubdomains[s].FcStar;
            }

            // Gather them in master
            //TODO: Perhaps I should cache them and reuse the unchanged ones.
            var transferrer = new TransferrerAltogetherFlattened(procs);
            Dictionary<int, Vector> allVectors_master = transferrer.GatherFromAllSubdomains(processVectors);

            // Calculate globalFcStar in master
            if (procs.IsMasterProcess)
            {
                var subdomainToVectors = allVectors_master.ChangeKey(model); //TODO: Changing the key of the Dictionary should be avoided

                // Give them to the global matrix manager so that it can create the global FcStar
                matrixManagerGlobal_master.CalcCoarseProblemRhs(subdomainToVectors); 
            }
        }

        public void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection)
        {
            // Calculate the subdomain KccStar in each process
            var processMatrices = new Dictionary<int, IMatrixView>();
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s] -> delegated to the SubdomainMatrixManager
                // globalKccStar = sum_over_s(Lc[s]^T * KccStar[s] * Lc[s]) -> delegated to the GlobalMatrixManager 
                ISubdomain subdomain = model.GetSubdomain(s);
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader + "Calculating Schur complement of remainder dofs"
                        + $" for the stiffness of subdomain {subdomain.ID}");
                    matrixManagerSubdomains[s].CalcCoarseProblemSubmatrices(); //TODO: At this point Kcc and Krc can be cleared. Maybe Krr too.
                }
                processMatrices[s] = matrixManagerSubdomains[s].CoarseProblemSubmatrix;
            }

            // Gather them in master
            //TODO: Perhaps I should cache them and reuse the unchanged ones.
            var transferrer = new TransferrerAltogetherFlattened(procs);
            Dictionary<int, IMatrixView> allMatrices_master = transferrer.GatherFromAllSubdomains(processMatrices);

            // Calculate globalKccStar and invert it in master
            if (procs.IsMasterProcess)
            {
                var subdomainsToMatrices = allMatrices_master.ChangeKey(model); //TODO: Changing the key of the Dictionary should be avoided

                // Give them to the global matrix manager so that it can create the global KccStar
                matrixManagerGlobal_master.CalcInverseCoarseProblemMatrix(cornerNodeSelection, subdomainsToMatrices);
            }
        }

        public void ClearCoarseProblemRhs()
        {
            if (procs.IsMasterProcess) matrixManagerGlobal_master.ClearCoarseProblemRhs();
        }

        public void ClearInverseCoarseProblemMatrix()
        {
            if (procs.IsMasterProcess) matrixManagerGlobal_master.ClearInverseCoarseProblemMatrix();
        }

        IFetiSubdomainMatrixManager IFetiMatrixManager.GetSubdomainMatrixManager(ISubdomain subdomain) 
            => GetFetiDPSubdomainMatrixManager(subdomain);

        public IFetiDPSubdomainMatrixManager GetFetiDPSubdomainMatrixManager(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return matrixManagerSubdomains[subdomain.ID];
        }

        public Vector MultiplyInverseCoarseProblemMatrix(Vector vector)
        {
            procs.CheckProcessIsMaster();
            return matrixManagerGlobal_master.MultiplyInverseCoarseProblemMatrixTimes(vector);
        }

        public DofPermutation ReorderGlobalCornerDofs()
        {
            procs.CheckProcessIsMaster();
            return matrixManagerGlobal_master.ReorderGlobalCornerDofs();
        }

        public DofPermutation ReorderSubdomainInternalDofs(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return matrixManagerSubdomains[subdomain.ID].ReorderInternalDofs();
        }

        public DofPermutation ReorderSubdomainRemainderDofs(ISubdomain subdomain)
        {
            procs.CheckProcessMatchesSubdomain(subdomain.ID);
            return matrixManagerSubdomains[subdomain.ID].ReorderRemainderDofs();
        }

        private Vector[] GatherCondensedRhsVectors(Vector subdomainVector)
        {
            //TODO: Perhaps I should cache them and reuse the unchanged ones. Use dedicated communication classes for this.
            return procs.Communicator.Gather(subdomainVector, procs.MasterProcess);
        }

        private IMatrixView[] GatherSchurComplementsOfRemainderDofs(IMatrixView subdomainMatrix)
        {
            //TODO: Perhaps I should cache them and reuse the unchanged ones. Use dedicated communication classes for this.
            return procs.Communicator.Gather(subdomainMatrix, procs.MasterProcess);
        }
    }
}
