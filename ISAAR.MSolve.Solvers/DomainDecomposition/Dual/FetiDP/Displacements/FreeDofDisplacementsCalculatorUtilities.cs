using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements
{
    internal static class FreeDofDisplacementsCalculatorUtilities
    {
        internal static void CalcAndStoreFreeDisplacements(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator, 
            IFetiDPMatrixManager matrixManager, ILagrangeMultipliersEnumerator lagrangesEnumerator, Vector lagranges, Vector uc)
        {
            IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            SignedBooleanMatrixColMajor Br = lagrangesEnumerator.GetBooleanMatrix(subdomain);

            // ur[s] = inv(Krr[s]) * (fr[s] - Br[s]^T * lagranges - Krc[s] * Bc[s] * uc)
            Vector BrLambda = Br.Multiply(lagranges, true);
            Vector KrcBcUc = Bc.Multiply(uc);
            KrcBcUc = subdomainMatrices.MultiplyKrcTimes(KrcBcUc);
            Vector ur = subdomainMatrices.Fr.Copy();
            ur.SubtractIntoThis(BrLambda);
            ur.SubtractIntoThis(KrcBcUc);
            ur = subdomainMatrices.MultiplyInverseKrrTimes(ur);

            // uf[s] = union(ur[s], ubc[s])
            // Remainder dofs
            var uf = Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs);
            int[] remainderDofs = dofSeparator.GetRemainderDofIndices(subdomain);
            uf.CopyNonContiguouslyFrom(remainderDofs, ur);

            // Corner dofs: ubc[s] = Bc[s] * uc
            Vector ubc = Bc.Multiply(uc);
            int[] cornerDofs = dofSeparator.GetCornerDofIndices(subdomain);
            uf.CopyNonContiguouslyFrom(cornerDofs, ubc);

            // Store uf[s]
            subdomainMatrices.LinearSystem.SolutionConcrete = uf;
        }
    }
}
