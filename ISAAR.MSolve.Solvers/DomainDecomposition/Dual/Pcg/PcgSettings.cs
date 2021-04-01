using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg
{
    public class PcgSettings
    {
        public IFetiPcgConvergenceFactory ConvergenceStrategyFactory { get; set; } =
            new ApproximateResidualConvergence.Factory();
        public double ConvergenceTolerance { get; set; } = 1E-7;
        public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);
    }
}
