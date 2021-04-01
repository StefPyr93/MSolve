using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    /// <summary>
    /// Subdoomain state update class that accounts for non zero initial conditions (displacements).
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class NonLinearSubdomainUpdaterWithInitialConditions : INonLinearSubdomainUpdater
	{
		private readonly Subdomain subdomain;

		public NonLinearSubdomainUpdaterWithInitialConditions(Subdomain subdomain)
		{
			this.subdomain = subdomain;
		}

		public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, Node> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements, ref BooleanArray isNodeUpdated, IVectorView subdomainLinearSystemSolution, ref BooleanArray areBoundaryNodesUpdated) //TODO leave 
		{
			if (!CnstValues.useV2FiniteElements)
			{
				return this.subdomain.GetRHSFromSolutionWithInitialDisplacemntsEffect(solution, dSolution, boundaryNodes,
				 initialConvergedBoundaryDisplacements, totalBoundaryDisplacements,
				 nIncrement, totalIncrements);
			}
			else
			{
				return this.subdomain.GetRHSFromSolutionWithInitialDisplacemntsEffect(solution, dSolution, boundaryNodes,
				 initialConvergedBoundaryDisplacements, totalBoundaryDisplacements,
				 nIncrement, totalIncrements, ref isNodeUpdated, subdomainLinearSystemSolution, ref areBoundaryNodesUpdated);
			}
		}



		public void ResetState()
		{
			this.subdomain.ClearMaterialStresses();
		}

		public void UpdateState()
		{
			this.subdomain.SaveMaterialState();
		}

		public void ScaleConstraints(double scalingFactor)
		{
			throw new NotSupportedException();
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution) //TODO leave 
		{
			throw new NotSupportedException();
			return this.subdomain.GetRhsFromSolution(solution, dSolution);
		}
	}
}
