using System;
using System.Collections.Generic;
using System.Text;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;

namespace MGroup.Stochastic.Interfaces
{
    public interface IProbabilityDistribution<Tobservations> : IMultivariateDistribution<Tobservations>, ISampleableDistribution<Tobservations>
    {
    }
}
