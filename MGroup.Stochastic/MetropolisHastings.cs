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
    public class MetropolisHastings : MetropolisHastings<double, MultivariateNormalDistribution>
    {


        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings"/> algorithm.
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="target">A function specifying the PDF or the log-PDF of the distribution to be sampled.</param>
        /// <param name="proposal">The proposal distribution that is used to generate new parameter samples to be explored.</param>
        /// 
        public MetropolisHastings(int dimensions, Func<double[], double> target, MultivariateNormalDistribution proposal)
            : base(dimensions, target, proposal)
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings"/> algorithm.
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="target">A function specifying the PDF or the log-PDF of the distribution to be sampled.</param>
        /// 
        public MetropolisHastings(int dimensions, Func<double[], double> target)
            : base(dimensions, target, new MultivariateNormalDistribution(dimensions))
        {
        }


        public MetropolisHastings(int dimensions, Func<double[], double> prior, Func<double[], double> likelihood, MultivariateNormalDistribution proposal)
            : base(dimensions, prior, likelihood, proposal)  //new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()))
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings"/> algorithm.
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="prior">A function specifying the PDF or the log-PDF of the prior distribution.</param>
        /// <param name="likelihood">A function specifying the PDF or the log-PDF of the likelihood function.</param>
        /// 
        public MetropolisHastings(int dimensions, Func<double[], double> prior, Func<double[], double> likelihood)
            : base(dimensions, prior, likelihood, new MultivariateNormalDistribution(dimensions))  //new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()))
        {
        }

        /// <summary>
        ///   Creates a new <see cref="MetropolisHastings"/> sampler using independent Normal distributions
        ///   as the parameter proposal generation priors. 
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="target">A function specifying the PDF or the log-PDF of the distribution to be sampled.</param>
        /// 
        /// <returns>A sampling algorithm that can generate samples from the target distribution.</returns>
        /// 
        public static MetropolisHastings<double, Independent<NormalDistribution>> Continuous(int dimensions, Func<double[], double> target)
        {
            return new MetropolisHastings<double, Independent<NormalDistribution>>(dimensions, target,
                new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()));
        }

        /// <summary>
        ///   Creates a new <see cref="MetropolisHastings"/> sampler using independent Normal distributions
        ///   as the parameter proposal generation priors. 
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="distribution">The target distribution whose samples should be generated.</param>
        /// 
        /// <returns>A sampling algorithm that can generate samples from the target distribution.</returns>
        /// 
        public static MetropolisHastings<double, Independent<NormalDistribution>, T> Continuous<T>(int dimensions, T distribution)
            where T : IMultivariateDistribution<double[]>
        {
            return new MetropolisHastings<double, Independent<NormalDistribution>, T>(distribution,
                new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()));
        }

        /// <summary>
        ///   Creates a new <see cref="MetropolisHastings"/> sampler using symmetric geometric distributions
        ///   as the parameter proposal generation priors. 
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="target">A function specifying the PDF or the log-PDF of the distribution to be sampled.</param>
        /// 
        /// <returns>A sampling algorithm that can generate samples from the target distribution.</returns>
        /// 
        public static MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>> Discrete(int dimensions, Func<int[], double> target)
        {
            return new MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>>(dimensions, target,
                new Independent<SymmetricGeometricDistribution, int>(dimensions, () => new SymmetricGeometricDistribution(0.5)));
        }

        /// <summary>
        ///   Creates a new <see cref="MetropolisHastings"/> sampler using symmetric geometric distributions
        ///   as the parameter proposal generation priors. 
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="distribution">The target distribution whose samples should be generated.</param>
        /// 
        /// <returns>A sampling algorithm that can generate samples from the target distribution.</returns>
        /// 
        public static MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>, T> Discrete<T>(int dimensions, T distribution)
            where T : IMultivariateDistribution<int[]>
        {
            return new MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>, T>(distribution,
                new Independent<SymmetricGeometricDistribution, int>(dimensions, () => new SymmetricGeometricDistribution(0.5)));
        }

    }

}
