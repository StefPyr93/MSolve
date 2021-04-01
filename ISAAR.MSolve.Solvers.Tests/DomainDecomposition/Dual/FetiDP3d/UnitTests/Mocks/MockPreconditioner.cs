using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks
{
    public class MockPreconditioner : IFetiPreconditioner
    {
        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.CopyFrom(Example4x4x4Quads.ExpectedGlobalMatrices.PreconditionerLumped * rhs);
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs)
        {
            lhs.CopyFrom(Example4x4x4Quads.ExpectedGlobalMatrices.PreconditionerLumped * rhs);
        }
    }
}
