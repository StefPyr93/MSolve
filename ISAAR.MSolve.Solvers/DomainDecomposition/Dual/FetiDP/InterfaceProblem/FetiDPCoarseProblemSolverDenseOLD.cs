using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public class FetiDPCoarseProblemSolverDenseOLD : IFetiDPCoarseProblemSolverOLD
    {
        private readonly IModel model;
        private Matrix inverseGlobalKccStar;

        public FetiDPCoarseProblemSolverDenseOLD(IModel model)
        {
            this.model = model;
        }

        public void ClearCoarseProblemMatrix()
        {
            inverseGlobalKccStar = null;
        }

        public void CreateAndInvertCoarseProblemMatrix(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains,
            FetiDPDofSeparatorOLD dofSeparator, Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers)
        {
            this.inverseGlobalKccStar = CreateGlobalKccStar(dofSeparator, matrixManagers);
            inverseGlobalKccStar.InvertInPlace();
        }

        public Vector CreateCoarseProblemRhs(FetiDPDofSeparatorOLD dofSeparator, 
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers, 
            Dictionary<int, Vector> fr, Dictionary<int, Vector> fbc)
        {
            // Static condensation for the force vectors
            var globalFcStar = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs);
            for (int s = 0; s < model.NumSubdomains; ++s)
            {
                IFetiDPSubdomainMatrixManagerOLD matrices = matrixManagers[s];

                // fcStar[s] = fbc[s] - Krc[s]^T * inv(Krr[s]) * fr[s]
                // globalFcStar = sum_over_s(Lc[s]^T * fcStar[s])
                UnsignedBooleanMatrix Lc = dofSeparator.CornerBooleanMatrices[s];
                Vector temp = matrices.MultiplyInverseKrrTimes(fr[s]);
                temp = matrices.MultiplyKcrTimes(temp);
                Vector fcStar = fbc[s] - temp;
                globalFcStar.AddIntoThis(Lc.Multiply(fcStar, true));
            }
            return globalFcStar;
        }

        public Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector) => inverseGlobalKccStar * vector;

        public void ReorderCornerDofs(FetiDPDofSeparatorOLD dofSeparator)
        {
            // Do nothing, since the sparsity pattern is irrelevant for dense matrices.
        }

        private Matrix CreateGlobalKccStar(FetiDPDofSeparatorOLD dofSeparator,
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers)
        {
            // Static condensation of remainder dofs (Schur complement).
            var globalKccStar = Matrix.CreateZero(dofSeparator.NumGlobalCornerDofs, dofSeparator.NumGlobalCornerDofs);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                int s = subdomain.ID;
                IFetiDPSubdomainMatrixManagerOLD matrices = matrixManagers[s];

                // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s]
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Calculating Schur complement of remainder dofs"
                        + " for the stiffness of subdomain {s}");
                    matrices.CalcSchurComplementOfRemainderDofs(); //TODO: At this point Kcc and Krc can be cleared. Maybe Krr too.
                }

                // globalKccStar = sum_over_s(Lc[s]^T * KccStar[s] * Lc[s])
                UnsignedBooleanMatrix Lc = dofSeparator.CornerBooleanMatrices[s];
                globalKccStar.AddIntoThis(Lc.ThisTransposeTimesOtherTimesThis(matrices.SchurComplementOfRemainderDofs));
            }
            return globalKccStar;
        }
    }
}
