using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements
{
    public class FreeDofDisplacementsCalculatorSerial : IFreeDofDisplacementsCalculator
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiDPMatrixManager matrixManager;
        private readonly IModel model;

        public FreeDofDisplacementsCalculatorSerial(IModel model, IFetiDPDofSeparator dofSeparator,
            IFetiDPMatrixManager matrixManager, ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.matrixManager = matrixManager;
            this.lagrangesEnumerator = lagrangesEnumerator;
        }

        public void CalculateSubdomainDisplacements(Vector lagranges, IFetiDPFlexibilityMatrix flexibility)
        {
            Vector uc = CalcCornerDisplacements(flexibility, lagranges);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                FreeDofDisplacementsCalculatorUtilities.CalcAndStoreFreeDisplacements(subdomain, dofSeparator, matrixManager,
                    lagrangesEnumerator, lagranges, uc);
            }
        }

        private Vector CalcCornerDisplacements(IFetiDPFlexibilityMatrix flexibility, Vector lagranges)
        {
            // uc = inv(KccStar) * (fcStar + FIrc^T * lagranges)
            Vector temp = flexibility.MultiplyFIrcTransposed(lagranges);
            temp.AddIntoThis(matrixManager.CoarseProblemRhs);
            return matrixManager.MultiplyInverseCoarseProblemMatrix(temp);
        }
    }
}
