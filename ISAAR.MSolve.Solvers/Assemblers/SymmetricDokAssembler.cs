using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;

//TODO: Instead of storing the raw CSC arrays, use a reusable DOK or SymmCscIndexer class. That class should provide methods to 
//      assemble the values part of the global matrix more efficiently than the general purpose DOK. The general purpose DOK 
//      should only be used to assemble the first global matrix and whenever the dof ordering changes. Now it is used everytime 
//      and the indexing arrays are discarded.
//TODO: I could also cache the symbolic factorization of SuiteSparse and reuse it. That would really speed up things.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in symmetric DOK format, namely only 
    /// the upper triangle is explicitly stored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricDokAssembler /*: IGlobalMatrixAssembler<SymmetricCscMatrix>*/
    {
        private const string name = "SymmetricDokAssembler"; // for error messages
        private ConstrainedMatricesAssembler constrainedAssembler = new ConstrainedMatricesAssembler();

        public SymmetricDokAssembler()
        {
        }

        public DokSymmetric BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var subdomainMatrix = DokSymmetric.CreateEmpty(numFreeDofs);

            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofIndices, subdomainDofIndices);
            }

            return subdomainMatrix;
        }

        public (DokSymmetric matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = freeDofOrdering.NumFreeDofs;
            var subdomainMatrix = DokSymmetric.CreateEmpty(numFreeDofs);

            //TODO: also reuse the indexers of the constrained matrices.
            constrainedAssembler.InitializeNewMatrices(freeDofOrdering.NumFreeDofs, constrainedDofOrdering.NumConstrainedDofs);

            // Process the stiffness of each element
            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofsFree, int[] subdomainDofsFree) = freeDofOrdering.MapFreeDofsElementToSubdomain(element);
                (int[] elementDofsConstrained, int[] subdomainDofsConstrained) =
                    constrainedDofOrdering.MapConstrainedDofsElementToSubdomain(element);

                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofsFree, subdomainDofsFree);
                constrainedAssembler.AddElementMatrix(elementMatrix, elementDofsFree, subdomainDofsFree,
                    elementDofsConstrained, subdomainDofsConstrained);
            }

            // Create the free and constrained matrices. 
            (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) = constrainedAssembler.BuildMatrices();
            return (subdomainMatrix, matrixConstrFree.TransposeToCSC(false), matrixConstrFree, matrixConstrConstr);
        }

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: Implement this
        }
    }
}
