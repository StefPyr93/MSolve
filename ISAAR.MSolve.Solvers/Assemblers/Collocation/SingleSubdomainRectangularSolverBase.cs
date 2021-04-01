using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public abstract class SingleSubdomainRectangularSolverBase<TMatrix>:ISolver
	where TMatrix: class,IMatrix
	{
        protected readonly IGlobalMatrixRectangularAssembler<TMatrix> assembler;
        protected readonly IStructuralAsymmetricModel model;
        protected readonly string name;
        protected readonly IAsymmetricSubdomain subdomain;
        protected readonly SingleSubdomainSystem<TMatrix> linearSystem;
        protected readonly IAsymmetricDofOrderer dofRowOrderer;
        protected readonly IDofOrderer dofColOrderer;

		protected SingleSubdomainRectangularSolverBase(IStructuralAsymmetricModel model, IAsymmetricDofOrderer dofRowOrderer,
			IDofOrderer dofColOrderer, IGlobalMatrixRectangularAssembler<TMatrix> assembler, string name)
		{
			if (model.NumSubdomains != 1) throw new InvalidSolverException(
				$"{name} can be used if there is only 1 subdomain");
			this.model = model;
			subdomain= model.Subdomains[0];
			linearSystem= new SingleSubdomainSystem<TMatrix>(subdomain);
			LinearSystems= new Dictionary<int, ILinearSystem>(){{subdomain.ID, linearSystem}};
			linearSystem.MatrixObservers.Add(this);

			this.dofRowOrderer = dofRowOrderer;
			this.dofColOrderer = dofColOrderer;
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
            TMatrix matrix = assembler.BuildGlobalMatrix(subdomain.FreeDofRowOrdering, subdomain.FreeDofColOrdering,
                    subdomain.EnumerateElements(), elementMatrixProvider);
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);
            return new Dictionary<int, IMatrix> { { subdomain.ID, matrix } };
        }
        

        public virtual Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView
            matrixConstrFree,
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
                subdomain.FreeDofRowOrdering, subdomain.FreeDofColOrdering, subdomain.ConstrainedDofRowOrdering,
                subdomain.ConstrainedDofColOrdering, subdomain.EnumerateElements(), elementMatrixProvider);
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);
            return new Dictionary<int, (IMatrix, IMatrixView, IMatrixView, IMatrixView)>
            {
                { subdomain.ID, (Aff, Afc, Acf, Acc) }
            };
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            if (otherMatrix.Count!=1) throw new InvalidSolverException("There can only be 1 subdomain when using this solver");
            KeyValuePair<int, IMatrixView> idMatrixPair = otherMatrix.First();
            int id = idMatrixPair.Key;
            Debug.Assert(id == subdomain.ID,
                "The matrix that will be multiplied with the inverse system matrix belongs to a different subdomain.");
            Matrix result = InverseSystemMatrixTimesOtherMatrix(idMatrixPair.Value);
            return new Dictionary<int, Matrix>{{id, result}};
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
		{
            var watch = new Stopwatch();
            watch.Start();

            assembler.HandleDofOrderingWillBeModified();

            //TODO xwsimo: Update AsymmetricDofOrderer so that it does the following itself, like DofOrderer does. 
            IGlobalFreeDofOrdering previousOrdering = model.GlobalDofOrdering; // will be overwritten
            dofColOrderer.OrderFreeDofs(model);
            IGlobalFreeDofOrdering globalRowOrdering = dofRowOrderer.OrderFreeDofs(model);

            model.GlobalRowDofOrdering = globalRowOrdering;
            model.GlobalColDofOrdering = model.GlobalDofOrdering;

            foreach (var subdomain in model.Subdomains)
            {
                subdomain.FreeDofRowOrdering = globalRowOrdering.GetSubdomainDofOrdering(subdomain);
                subdomain.FreeDofColOrdering = model.GlobalColDofOrdering.GetSubdomainDofOrdering(subdomain);

                if (alsoOrderConstrainedDofs)
                {
                    subdomain.ConstrainedDofRowOrdering = dofRowOrderer.OrderConstrainedDofs(subdomain);
                    subdomain.ConstrainedDofColOrdering = dofColOrderer.OrderConstrainedDofs(subdomain);
                }
            }

            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            Logger.LogNumDofs("Global rows", globalRowOrdering.NumGlobalFreeDofs);
            Logger.LogNumDofs("Global columns", model.GlobalColDofOrdering.NumGlobalFreeDofs);
        }

        public abstract void Initialize();
        public abstract void HandleMatrixWillBeSet();
        public abstract void PreventFromOverwrittingSystemMatrices();
        public abstract void Solve();
        protected abstract Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix);
    }
}
