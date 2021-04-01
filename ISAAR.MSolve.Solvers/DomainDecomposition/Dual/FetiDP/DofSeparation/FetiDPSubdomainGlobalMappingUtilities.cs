using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    internal static class FetiDPSubdomainGlobalMappingUtilities
    {
        internal static void ScaleSubdomainFreeDisplacements(IFetiDPDofSeparator dofSeparator, 
            IStiffnessDistribution distribution, ISubdomain subdomain, Vector subdomainFreeDisplacements)
        {
            // Boundary remainder dofs: Scale them so that they can be just added at global level
            int[] remainderToSubdomainDofs = dofSeparator.GetRemainderDofIndices(subdomain);
            double[] boundaryDofCoeffs = distribution.CalcBoundaryDofCoefficients(subdomain);
            int[] boundaryDofIndices = dofSeparator.GetBoundaryDofIndices(subdomain);
            for (int i = 0; i < boundaryDofIndices.Length; ++i)
            {
                int idx = remainderToSubdomainDofs[boundaryDofIndices[i]];
                subdomainFreeDisplacements[idx] *= boundaryDofCoeffs[i];
            }

            // Boundary corner dofs: Scale them based on their multiplicity so that they can be just added at global level
            //TODO: Is the multiplicity scaling for corner dofs correct? I also used it when distributing the external loads.
            int[] cornerToSubdomainDofs = dofSeparator.GetCornerDofIndices(subdomain);
            IReadOnlyList<(INode node, IDofType dofType)> cornerDofs = dofSeparator.GetCornerDofs(subdomain);
            for (int i = 0; i < cornerToSubdomainDofs.Length; ++i)
            {
                int idx = cornerToSubdomainDofs[i];
                subdomainFreeDisplacements[idx] /= cornerDofs[i].node.Multiplicity;
            }
        }
    }
}
