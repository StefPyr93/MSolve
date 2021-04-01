using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: not sure this class should be in this namespace
namespace ISAAR.MSolve.Solvers
{
    /// <summary>
    /// Base implementation for solver that do not use domain decomposition.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix">The type of the linear system's matrix.</typeparam>
    public abstract class SingleSubdomainSolverBase<TMatrix> : ISolver
        where TMatrix : class, IMatrix
    {
        protected readonly IGlobalMatrixAssembler<TMatrix> assembler;
        protected readonly IDofOrderer dofOrderer;
        protected readonly IModel model;
        protected readonly string name; // for error messages
        protected readonly ISubdomain subdomain;
        protected readonly SingleSubdomainSystem<TMatrix> linearSystem;

        protected SingleSubdomainSolverBase(IModel model, IDofOrderer dofOrderer, 
            IGlobalMatrixAssembler<TMatrix> assembler, string name)
        {
            if (model.NumSubdomains != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            subdomain = model.EnumerateSubdomains().First();

            linearSystem = new SingleSubdomainSystem<TMatrix>(subdomain);
            LinearSystems = new Dictionary<int, ILinearSystem>() { { subdomain.ID, linearSystem } };
            linearSystem.MatrixObservers.Add(this);

            this.dofOrderer = dofOrderer;
            this.assembler = assembler;
            this.Logger = new SolverLoggerOLD(name);
        }

        public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }
        public SolverLoggerOLD Logger { get; }
        public string Name { get; }

        //TODO: This should be defined only for DDM solvers
        public INodalLoadDistributor NodalLoadDistributor => throw new NotImplementedException();

        public virtual Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var watch = new Stopwatch();
            watch.Start();
            TMatrix matrix = assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, 
                subdomain.EnumerateElements(), elementMatrixProvider);
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);
            return new Dictionary<int, IMatrix> { { subdomain.ID, matrix } };
        }

        public virtual Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            var watch = new Stopwatch();
            watch.Start();
            if (subdomain.ConstrainedDofOrdering == null)
            {
                throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs,"
                    + " they must have been ordered first.");
            }
            (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc) = assembler.BuildGlobalSubmatrices(
                subdomain.FreeDofOrdering, subdomain.ConstrainedDofOrdering, subdomain.EnumerateElements(), elementMatrixProvider);
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);
            return new Dictionary<int, (IMatrix, IMatrixView, IMatrixView, IMatrixView)>
            {
                { subdomain.ID, (Aff, Afc, Acf, Acc) }
            };
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            if (otherMatrix.Count != 1) throw new InvalidSolverException("There can only be 1 subdomain when using this solver");
            KeyValuePair<int, IMatrixView> idMatrixPair = otherMatrix.First();
            int id = idMatrixPair.Key;
            Debug.Assert(id == subdomain.ID, 
                "The matrix that will be multiplied with the inverse system matrix belongs to a different subdomain.");
            Matrix result = InverseSystemMatrixTimesOtherMatrix(idMatrixPair.Value);
            return new Dictionary<int, Matrix>() { { id, result } };
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
        {
            var watch = new Stopwatch();
            watch.Start();

            assembler.HandleDofOrderingWillBeModified();
            dofOrderer.OrderFreeDofs(model);

            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();

            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            Logger.LogNumDofs("Global dofs", model.GlobalDofOrdering.NumGlobalFreeDofs);
        }

        public abstract void Initialize();
        public abstract void HandleMatrixWillBeSet();
        public abstract void PreventFromOverwrittingSystemMatrices();
        public abstract void Solve();
        protected abstract Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix);
    }
}
