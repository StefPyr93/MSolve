using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: What about transpose(this) * other * this operations? More importantly, what about accessing rows2cols (or cols2rows) map? 
//TODO: Similarly, what about transpose(this) * factorizedOther * this operations?
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    public interface IMappingMatrix : IBounded2D
    {
        Matrix CopyToFullMatrix(); //TODO: Should this be here?

        Vector Multiply(Vector vector, bool transposeThis = false);

        Matrix MultiplyRight(Matrix other, bool transposeThis = false);
    }
}
