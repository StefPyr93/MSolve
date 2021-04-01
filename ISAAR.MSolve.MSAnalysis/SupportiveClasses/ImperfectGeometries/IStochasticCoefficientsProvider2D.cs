using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticCoefficientsProvider2D //: IStochasticCoefficientsProvider
    {
        double[] RandomVariables { get; set; }
        double GetCoefficient(double meanValue, double[] coordinates);
        double[] GetDerivative(double[] coordinates);
    }
}
