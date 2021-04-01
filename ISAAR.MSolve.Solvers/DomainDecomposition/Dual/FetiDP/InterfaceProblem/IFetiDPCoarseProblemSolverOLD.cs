using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPCoarseProblemSolverOLD
    {
        void ClearCoarseProblemMatrix();

        Vector CreateCoarseProblemRhs(FetiDPDofSeparatorOLD dofSeparator,
            Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers,
            Dictionary<int, Vector> fr, Dictionary<int, Vector> fbc);

        //TODO: Perhaps corner nodes of each subdomain should be stored in FetiDPDofSeparator.
        void CreateAndInvertCoarseProblemMatrix(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains, 
            FetiDPDofSeparatorOLD dofSeparator, Dictionary<int, IFetiDPSubdomainMatrixManagerOLD> matrixManagers);

        Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector);

        void ReorderCornerDofs(FetiDPDofSeparatorOLD dofSeparator);
    }
}
