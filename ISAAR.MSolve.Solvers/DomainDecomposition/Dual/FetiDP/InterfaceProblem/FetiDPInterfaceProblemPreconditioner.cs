using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

//TODO: Perhaps remove this class and make IFetiPreconditioner implement IPreconditioner.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public class FetiDPInterfaceProblemPreconditioner : IPreconditioner
    {
        private readonly IFetiPreconditioner fetiPreconditioner;

        internal FetiDPInterfaceProblemPreconditioner(IFetiPreconditioner fetiPreconditioner)
        {
            this.fetiPreconditioner = fetiPreconditioner;
        }

        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
            //      IVectorView and IVector
            var lhs = (Vector)lhsVector;//bookmark1
            var rhs = (Vector)rhsVector;
            fetiPreconditioner.SolveLinearSystem(rhs, lhs);
            if (CnstValues.printPreconditoner)
            {
                var prec =Matrix.CreateZero(lhs.Length, rhs.Length);
                fetiPreconditioner.SolveLinearSystems(Matrix.CreateIdentity(rhs.Length),prec);
                string pathprec = (new CnstValues()).solverPath + @"\preconditioner.txt";
                (new FullMatrixWriter()).WriteToFile(prec, pathprec);
                CnstValues.printPreconditoner = false;
            }
            if (CnstValues.printNRiterPreconditioner == CnstValues.analyzerNRIter )
            {
                var prec = Matrix.CreateZero(lhs.Length, rhs.Length);
                fetiPreconditioner.SolveLinearSystems(Matrix.CreateIdentity(rhs.Length), prec);
                string pathprec = (new CnstValues()).solverPath + @"\preconditioner.txt";
                (new FullMatrixWriter()).WriteToFile(prec, pathprec);
                CnstValues.printPreconditoner = false;
            }
        }
    }
}
