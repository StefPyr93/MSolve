﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Analysis
{
    class LinearStaticAnalysisCSC: ILinearStaticAnalysis
    {
        private readonly Model2D model;

        public LinearStaticAnalysisCSC(Model2D model)
        {
            this.model = model;
        }

        public IDOFEnumerator DOFEnumerator { get; private set; }
        public Vector Solution { get; private set; }

        public void Solve()
        {
            DOFEnumerator = new DOFEnumerator(model);
            (DOKSymmetricColMajor matrix, Vector rhs) = ReduceToSimpleLinearSystem();
            using (CholeskySuiteSparse factorization = matrix.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                Solution = factorization.SolveLinearSystem(rhs);
            }
        }

        public void PrintSolution()
        {
            Console.WriteLine("Displacements: ");
            for (int n = 0; n < model.Nodes.Count; ++n)
            {
                int xDof = DOFEnumerator.GetFreeDofOf(model.Nodes[n], DisplacementDOF.X);
                double dx = (xDof < 0) ? 0 : Solution[xDof];
                int yDof = DOFEnumerator.GetFreeDofOf(model.Nodes[n], DisplacementDOF.Y);
                double dy = (yDof < 0) ? 0 : Solution[yDof];

                Console.WriteLine("Node " + n + ": dx = " + dx + "\t\t , dy = " + dy);
            }
        }

        private (DOKSymmetricColMajor matrix, Vector rhs) ReduceToSimpleLinearSystem()
        {
            /// The extended linear system is:
            /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
            /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
            /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
            /// To solve the system (for the unknowns ul):
            /// i) Kuu * uu = Fu - Kuc * uc = Feff
            /// ii) uu = Kuu \ Feff
            /// 
            (DOKSymmetricColMajor Kuu, Matrix Kuc) = SingleGlobalDOKAssembler.BuildGlobalMatrix(model, DOFEnumerator);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            Vector Fu = model.CalculateFreeForces(DOFEnumerator); //TODO: fix MKL dlls
            Vector uc = model.CalculateConstrainedDisplacements(DOFEnumerator);
            Vector Feff = Fu - Kuc * uc;
            return (Kuu, Feff);
        }
    }
}
