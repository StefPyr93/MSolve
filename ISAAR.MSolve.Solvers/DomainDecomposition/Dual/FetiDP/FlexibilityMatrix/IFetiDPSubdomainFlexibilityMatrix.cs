using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix
{
    public interface IFetiDPSubdomainFlexibilityMatrix
    {
        Vector MultiplyFIrc(Vector vector);

        Vector MultiplyFIrcTransposed(Vector vector);

        Vector MultiplyFIrr(Vector vector);

        (Vector FIrrTimesVector, Vector FIrcTransposedTimesVector) MultiplyFIrrAndFIrcTransposedTimesVector(Vector vector);
    }
}
