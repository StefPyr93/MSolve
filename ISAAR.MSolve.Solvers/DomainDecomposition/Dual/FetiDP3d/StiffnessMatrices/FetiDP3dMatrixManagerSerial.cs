using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Unify this with FetiDPMatrixManagerSerial. SubdomainMatrixManager should provide access to a list of matrices, not just
//      KccStar.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices
{
    public class FetiDP3dMatrixManagerSerial : IFetiDP3dMatrixManager
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiDPGlobalMatrixManager matrixManagerGlobal;
        private readonly Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager> matrixManagersSubdomain;
        private readonly IModel model;
        private readonly string msgHeader;

        public FetiDP3dMatrixManagerSerial(IModel model, IFetiDPDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IAugmentationConstraints augmentationConstraints, 
            IFetiDP3dMatrixManagerFactory matrixManagerFactory)
        {
            this.model = model;
            this.lagrangesEnumerator = lagrangesEnumerator;
            this.augmentationConstraints = augmentationConstraints;
            this.matrixManagersSubdomain = new Dictionary<ISubdomain, IFetiDPSubdomainMatrixManager>();
            matrixManagerGlobal = matrixManagerFactory.CreateGlobalMatrixManager(model, dofSeparator, augmentationConstraints);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                this.matrixManagersSubdomain[sub] = 
                    matrixManagerFactory.CreateSubdomainMatrixManager(sub, dofSeparator, lagrangesEnumerator,
                    augmentationConstraints);
            }
            this.msgHeader = $"{this.GetType().Name}: ";
        }

        //TODO: These 2 vectors should be stored and calculated by the class that manages global vectors, such as fcStar. 
        //      For now this component is IFetiDP3dGlobalMatrixManager.
        public Vector CoarseProblemRhs { get; private set; }
        public Vector GlobalDr { get; private set; }

        public void CalcCoarseProblemRhs()
        {
            // fcStarTilde = [fcStar ; -Qr^T * dr]
            CalcGlobalFcStar();
            CalcGlobalDr();
            Vector QrDr = augmentationConstraints.MatrixGlobalQr.Multiply(GlobalDr.Scale(-1), true);
            Vector fcStarTilde = matrixManagerGlobal.CoarseProblemRhs.Append(QrDr);
            this.CoarseProblemRhs = fcStarTilde;
        }

        public void CalcInverseCoarseProblemMatrix(ICornerNodeSelection cornerNodeSelection)
        {
            // Calculate KccStarTilde of each subdomain
            // globalKccStarTilde = sum_over_s(Lc[s]^T * KccStarTilde[s] * Lc[s]) -> delegated to the GlobalMatrixManager 
            // Here we will just prepare the data for GlobalMatrixManager
            var allKccStarTilde = new Dictionary<ISubdomain, IMatrixView>();
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                IFetiDPSubdomainMatrixManager matrices = matrixManagersSubdomain[sub];
                if (sub.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader + "Calculating Schur complement of remainder dofs"
                        + $" for the stiffness of subdomain {sub.ID}");
                    matrices.CalcCoarseProblemSubmatrices(); //TODO: At this point Kcc and Krc can be cleared. Maybe Krr too.
                }
                allKccStarTilde[sub] = matrices.CoarseProblemSubmatrix;
            }

            // Give them to the global matrix manager so that it can create the global KccStar
            matrixManagerGlobal.CalcInverseCoarseProblemMatrix(cornerNodeSelection, allKccStarTilde);
        }

        public void ClearCoarseProblemRhs() => matrixManagerGlobal.ClearCoarseProblemRhs();
        public void ClearInverseCoarseProblemMatrix() => matrixManagerGlobal.ClearInverseCoarseProblemMatrix();

        public IFetiSubdomainMatrixManager GetSubdomainMatrixManager(ISubdomain subdomain)
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

        //TODO: For now this is a side effect of CalcCoarseProblemRhs(). It should be called explicitly by whoever calls 
        //CalcCoarseProblemRhs().
        private void CalcGlobalDr()
        {
            GlobalDr = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                //TODO: This should be done by the same class as subdomain's fcStar
                Vector subdomainDr = FetiDPInterfaceProblemUtilities.CalcSubdomainDr(sub, this, lagrangesEnumerator);
                GlobalDr.AddIntoThis(subdomainDr); //TODO: This should be done by the same class as global fcStar
            }
        }

        private void CalcGlobalFcStar()
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
    }
}
