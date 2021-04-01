using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.LinearSystems
{
    public class SingleSubdomainSystemMpi<TMatrix> : ISingleSubdomainLinearSystemMpi
        where TMatrix : class, IIndexable2D
    {
        protected const int initialSize = int.MinValue;

        public SingleSubdomainSystemMpi(ISubdomain subdomain)
        {
            this.Subdomain = subdomain;
        }

        public TMatrix Matrix { get; set; }
        IIndexable2D ILinearSystemMpi.Matrix
        {
            get => Matrix;
            set
            {
                if ((value.NumRows != this.Size) || (value.NumColumns != this.Size))
                {
                    throw new NonMatchingDimensionsException("The provided matrix does not match the dimensions or pattern of"
                        + " this linear system. Make sure that it is initialization was delegated to this linear system"
                        + " after the latest dof ordering.");
                }
                foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();
                Matrix = (TMatrix)value;
            }
        }

        public HashSet<ISystemMatrixObserver> MatrixObservers { get; } = new HashSet<ISystemMatrixObserver>();

        public Vector RhsConcrete { get; set; }
        IVector ILinearSystemMpi.RhsVector
        {
            get => RhsConcrete;
            set
            {
                if (value.Length != this.Size)
                {
                    throw new NonMatchingDimensionsException("The provided vector does not match the dimensions or pattern of"
                        + " this linear system. Make sure that it is initialization was delegated to this linear system"
                        + " after the latest dof ordering.");
                }
                RhsConcrete = (Vector)value;
            }
        }

        public int Size { get; private set; } = initialSize;

        IVectorView ILinearSystemMpi.Solution => SolutionConcrete;
        public Vector SolutionConcrete { get; set; }

        public ISubdomain Subdomain { get; }

        IVector ILinearSystemMpi.CreateZeroVector() => CreateZeroVectorConcrete();
        public Vector CreateZeroVectorConcrete()
        {
            if (Size == initialSize) throw new InvalidOperationException(
                "The linear system size must be set before creating vectors. First of all order the subdomain freedom degrees.");
            return Vector.CreateZero(Size);
        }

        public virtual void Reset()
        {
            foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();

            // Override the method if memory needs to be disposed in a more complicated way.
            RhsConcrete = null;
            SolutionConcrete = null;
            Matrix = null;

            if (Subdomain.FreeDofOrdering == null)
            {
                throw new InvalidOperationException("The freedom degrees of a subdomain must"
                    + " be ordered before defining the size of its corresponding linear system.");
            }
            Size = Subdomain.FreeDofOrdering.NumFreeDofs;
        }
    }
}
