using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

//TODO: Qr matrix should be defined per subdomain. Crucial for MPI
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public interface IAugmentationConstraints
    {
        DofTable GlobalAugmentationDofOrdering { get; }

        /// <summary>
        /// Qr is a (nL x na) matrix where nL is the number of global lagrange multipliers and na is 
        /// <see cref="NumGlobalAugmentationConstraints"/>.
        /// </summary>
        IMappingMatrix MatrixGlobalQr { get; }

        IMidsideNodesSelection MidsideNodesSelection { get; }

        /// <summary>
        /// The number of extra constraints for the 3D problem. E.g. in "A scalable dual–primal domain decomposition method, 
        /// Farhat et al, 2000" it is proposed to add 3 constraints (X,Y,Z) at the middle of each boundary edge between 
        /// subdomains.
        /// </summary>
        int NumGlobalAugmentationConstraints { get; }

        void CalcAugmentationMappingMatrices();

        DofTable GetAugmentationDofOrdering(ISubdomain subdomain);


        /// <summary>
        /// The augmented constraints of the subdomain are in the same order as the nodes in 
        /// <see cref="IMidsideNodesSelection.GetMidsideNodesOfSubdomain(ISubdomain)"/> of <see cref="MidsideNodesSelection"/>.
        /// </summary>
        /// <param name="subdomain"></param>
        /// <returns></returns>
        GlobalToLocalBooleanMatrix GetMatrixBa(ISubdomain subdomain); //TODO: This should be an interface instead. Hard to define interface though.

        IMappingMatrix GetMatrixR1(ISubdomain subdomain);

        int GetNumAugmentationDofs(ISubdomain subdomain);
    }
}
