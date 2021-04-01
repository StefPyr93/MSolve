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
    public class MetropolisHastings<TObservation, TProposalDistribution>
        : MetropolisHastings<TObservation>
        where TProposalDistribution : ISampleableDistribution<TObservation[]>
    {
        private TProposalDistribution proposal;

        /// <summary>
        ///   Gets or sets the move proposal distribution.
        /// </summary>
        /// 
        public new TProposalDistribution Proposal
        {
            get { return proposal; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings"/> algorithm.
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="target">A function specifying the log PDF or the log-PDF of the distribution to be sampled.</param>
        /// <param name="proposal">The proposal distribution that is used to generate new parameter samples to be explored.</param>
        /// 
        public MetropolisHastings(int dimensions, Func<TObservation[], double> target, TProposalDistribution proposal)
        {
            Initialize(dimensions, target, proposal);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings"/> algorithm.
        /// </summary>
        /// 
        /// <param name="dimensions">The number of dimensions in each observation.</param>
        /// <param name="prior">A function specifying the PDF or the log-PDF of the prior distribution.</param>
        /// <param name="likelihood">A function specifying the PDF or the log-PDF of the likelihood function.</param>
        /// <param name="proposal">The proposal distribution that is used to generate new parameter samples to be explored.</param>
        /// 
        public MetropolisHastings(int dimensions, Func<TObservation[], double> prior, Func<TObservation[], double> likelihood, TProposalDistribution proposal)
        {
            Initialize(dimensions, prior, likelihood, proposal);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="MetropolisHastings{TObservation, TProposalDistribution}"/> class.
        /// </summary>
        /// 
        protected MetropolisHastings()
        {

        }

        /// <summary>
        ///   Initializes the algorithm.
        /// </summary>
        /// 
        protected void Initialize(int dimensions, Func<TObservation[], double> target, TProposalDistribution proposal)
        {
            this.proposal = proposal;
            Initialize(dimensions, target, generate);
        }

        protected void Initialize(int dimensions, Func<TObservation[], double> prior, Func<TObservation[], double> likelihood, TProposalDistribution proposal)
        {
            this.proposal = proposal;
            Initialize(dimensions, prior, likelihood, generate);
        }

        private TObservation[] generate(TObservation[] current, TObservation[] next)
        {
            return proposal.Generate(result: next);
        }

        public void UpdateProposalDistribution(TProposalDistribution proposal)
        {
            this.proposal = proposal;
        }

    }

}