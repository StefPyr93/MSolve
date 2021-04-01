using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements
{
    public interface IFreeDofDisplacementsCalculator
    {
        void CalculateSubdomainDisplacements(Vector lagranges, IFetiDPFlexibilityMatrix flexibility);
    }
}
