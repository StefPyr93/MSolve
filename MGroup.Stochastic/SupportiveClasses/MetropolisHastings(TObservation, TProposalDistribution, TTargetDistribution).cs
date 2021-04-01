using Accord.Math.Random;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using System;

namespace MGroup.Stochastic
{
    /// <summary>
    ///   Metropolis-Hasting sampling algorithm.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm">
    ///       Wikipedia, The Free Encyclopedia. Metropolis-Hastings algorithm. 
    ///       Available on: https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm </a></description></item>
    ///     <item><description><a href="http://en.wikipedia.org/wiki/Joint_probability_distribution">
    ///       Darren Wilkinson, Metropolis Hastings MCMC when the proposal and target have differing support.  
    ///       Available on: https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/ </a></description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    public class MetropolisHastings<TObservation, TProposalDistribution, TTargetDistribution>
        : MetropolisHastings<TObservation, TProposalDistribution>
        where TProposalDistribution : ISampleableDistribution<TObservation[]>
        where TTargetDistribution : IMultivariateDistribution<TObservation[]>
    {
        TTargetDistribution target;

        /// <summary>
        ///   Gets the target distribution whose samples must be generated.
        /// </summary>
        /// 
        public TTargetDistribution Target { get { return target; } }


        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings{TObservation, TProposalDistribution, TTargetDistribution}"/> class.
        /// </summary>
        /// 
        /// <param name="target">The target distribution whose samples should be generated.</param>
        /// <param name="proposal">The proposal distribution that is used to generate new parameter samples to be explored.</param>
        ///
        public MetropolisHastings(TTargetDistribution target, TProposalDistribution proposal)
        {
            this.target = target;
            this.Initialize(target.Dimension, target.LogProbabilityFunction, proposal);
        }
    }
}
