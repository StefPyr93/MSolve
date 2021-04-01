using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.LinearSystems
{
    public interface ILinearSystemMpi
    {
        /// <summary>
        /// The linear system matrix.
        /// </summary>
        IIndexable2D Matrix { get; set; }

        /// <summary>
        /// Objects that will be notified when <see cref="Matrix"/> is modified.
        /// </summary>
        HashSet<ISystemMatrixObserver> MatrixObservers { get; }

        /// <summary>
        /// The right hand side vector.
        /// </summary>
        IVector RhsVector { get; set; }

        /// <summary>
        /// The number of linear equations, which is equal to the number of unknowns and the dimensions of the matrix 
        /// and vectors.
        /// </summary>
        int Size { get; }

        /// <summary>
        /// The subdomain corresponding to this linear system.
        /// </summary>
        ISubdomain Subdomain { get; } //TODO: delete this once subdomains have been abstracted.

        /// <summary>
        /// The solution vector.
        /// </summary>
        IVectorView Solution { get; }

        /// <summary>
        /// Initializes a new vector with zero entries. Its pattern depends on the solver used. The freedom degrees must be 
        /// ordered before this method can be called. Special attention is needed if the freedom degrees change during 
        /// the analysis (e.g. adaptive FEM, XFEM). 
        /// </summary>
        IVector CreateZeroVector();

        /// <summary>
        /// Clears all data stored in this <see cref="ILinearSystem"/> instance and sets its <see cref="Size"/> equal to the
        /// number of freedom degrees of the corresponding <see cref="Subdomain"/>.
        /// </summary>
        void Reset();
    }
}
