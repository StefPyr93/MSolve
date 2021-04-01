using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices
{
    public class FetiDPMatrixManagerSerial : IFetiDPMatrixManager
    {
        private readonly IFetiDPGlobalMatrixManager matrixManagerGlobal;
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager> matrixManagersSubdomain;
        private readonly IModel model;
        private readonly string msgHeader;

        public FetiDPMatrixManagerSerial(IModel model, IFetiDPDofSeparator dofSeparator,
            IFetiDPMatrixManagerFactory matrixManagerFactory)
        {
            this.model = model;
            this.matrixManagersSubdomain = new Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager>();
            matrixManagerGlobal = matrixManagerFactory.CreateGlobalMatrixManager(model, dofSeparator);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                this.matrixManagersSubdomain[sub] = matrixManagerFactory.CreateSubdomainMatrixManager(sub, dofSeparator);
            }
            this.msgHeader = $"{this.GetType().Name}: ";
        }

        public Vector CoarseProblemRhs => matrixManagerGlobal.CoarseProblemRhs;

        public void CalcCoarseProblemRhs()
        {
            // Calculate FcStar of each subdomain
            // fcStar[s] = fbc[s] - Krc[s]^T * inv(Krr[s]) * fr[s] -> delegated to the SubdomainMatrixManager
            // globalFcStar = sum_over_s(Lc[s]^T * fcStar[s]) -> delegated to the GlobalMatrixManager
            var allFcStar = new Dictionary<ISubdomain, Vector>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                matrixManagersSubdomain[sub].CalcCoarseProblemRhsSubvectors();
                allFcStar[sub] = matrixManagersSubdomain[sub].FcStar;
            }

            // Give them to the global matrix manager so that it can create the global FcStar
            matrixManagerGlobal.CalcCoarseProblemRhs(allFcStar);
        }

        public void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection)
        {
            // Calculate KccStar of each subdomain
            // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s] -> delegated to the SubdomainMatrixManager
            // globalKccStar = sum_over_s(Lc[s]^T * KccStar[s] * Lc[s]) -> delegated to the GlobalMatrixManager 
            var allKccStar = new Dictionary<ISubdomain, IMatrixView>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                if (sub.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader + "Calculating Schur complement of remainder dofs"
                        + $" for the stiffness of subdomain {sub.ID}");
                    matrixManagersSubdomain[sub].CalcCoarseProblemSubmatrices(); //TODO: At this point Kcc and Krc can be cleared. Maybe Krr too.
                } 
                allKccStar[sub] = matrixManagersSubdomain[sub].CoarseProblemSubmatrix;
            }

            // Give them to the global matrix manager so that it can create the global KccStar
            matrixManagerGlobal.CalcInverseCoarseProblemMatrix(cornerNodeSelection, allKccStar);
        }

        public void ClearCoarseProblemRhs() => matrixManagerGlobal.ClearCoarseProblemRhs();
        public void ClearInverseCoarseProblemMatrix() => matrixManagerGlobal.ClearInverseCoarseProblemMatrix();

        IFetiSubdomainMatrixManager IFetiMatrixManager.GetSubdomainMatrixManager(ISubdomain subdomain)
            => matrixManagersSubdomain[subdomain];

        public IFetiDPSubdomainMatrixManager GetFetiDPSubdomainMatrixManager(ISubdomain subdomain) 
            => matrixManagersSubdomain[subdomain];

        public Vector MultiplyInverseCoarseProblemMatrix(Vector vector) 
            => matrixManagerGlobal.MultiplyInverseCoarseProblemMatrixTimes(vector);

        public DofPermutation ReorderGlobalCornerDofs() => matrixManagerGlobal.ReorderGlobalCornerDofs();

        public DofPermutation ReorderSubdomainInternalDofs(ISubdomain subdomain) 
            => matrixManagersSubdomain[subdomain].ReorderInternalDofs();

        public DofPermutation ReorderSubdomainRemainderDofs(ISubdomain subdomain)
            => matrixManagersSubdomain[subdomain].ReorderRemainderDofs();
    }
}
