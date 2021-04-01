using Accord.Math.Random;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using MGroup.Stochastic.Interfaces;
using System;

namespace MGroup.Stochastic
{

    public class MetropolisHastingsNew
    {

        Func<double[], double> targetPDF;
        Func<double[], double> priorPDF;
        Func<double[], double> likelihoodPDF;
        IProbabilityDistribution<double[]> proposal;

        double[] current;
        double[] candidate;

        double pCurrent;

        int dimensions;
        int discard = 0; // steps to discard
        long steps = 0; // steps so far
        long accepts = 0; // steps accepted

        bool initialized = false;

        public double[] ModelOutput { get; set; }
        public double[] CandidateSample { get; set; }

        public bool islogPDF = false;

        /// <summary>
        ///   Gets the last successfully generated observation.
        /// </summary>
        /// 
        public double[] CurrentSample { get { return current; } }

        /// <summary>
        ///   Gets or sets a factory method to create random number generators used in this instance.
        /// </summary>
        /// 
        public Random RandomSource { get; set; }

        /// <summary>
        ///   Gets the log-probability of the <see cref="CurrentSample">last successfully
        ///   generated sample</see>.
        /// </summary>
        /// 
        public double CurrentValue { get { return pCurrent; } }

        /// <summary>
        ///   Gets the log-probability density function of the target distribution.
        /// </summary>
        /// 
        public Func<double[], double> ProbabilityDensityFunction
        {
            get { return targetPDF; }
        }

        /// <summary>
        ///   Gets or sets the move proposal distribution.
        /// </summary>
        /// 
        public IProbabilityDistribution<double[]> Proposal
        {
            get { return proposal; }
        }

        /// <summary>
        ///   Gets the acceptance rate for the proposals generated
        ///   by the <see cref="Proposal">proposal distribution</see>.
        /// </summary>
        /// 
        public double AcceptanceRate
        {
            get { return this.accepts / (double)this.steps; }
        }

        /// <summary>
        ///   Gets the number of dimensions in each observation.
        /// </summary>
        /// 
        public int NumberOfInputs
        {
            get { return dimensions; }
        }

        /// <summary>
        ///   Gets or sets how many initial samples will get discarded as part
        ///   of the initial thermalization (warm-up, initialization) process.
        /// </summary>
        /// 
        public int Discard
        {
            get { return discard; }
            set { discard = value; }
        }

        public MetropolisHastingsNew(int dimensions, IMultivariateDistribution<double[]> target, IProbabilityDistribution<double[]> proposal)
        {
            targetPDF = target.LogProbabilityFunction;
            this.proposal = proposal;
            this.dimensions = dimensions;
            Initialize();
        }

        public MetropolisHastingsNew(int dimensions, IMultivariateDistribution<double[]> prior, IMultivariateDistribution<double[]> likelihood, IProbabilityDistribution<double[]> proposal)
        {
            priorPDF = prior.LogProbabilityFunction;
            likelihoodPDF = likelihood.LogProbabilityFunction;
            this.proposal = proposal;
            this.dimensions = dimensions;
            Initialize();
        }

        protected void Initialize()
        {
            this.current = new double[dimensions];
            this.candidate = new double[dimensions];
            if (likelihoodPDF == null)
                this.pCurrent = targetPDF(current);
            else
                this.pCurrent = priorPDF(current) * likelihoodPDF(ModelOutput);
            this.RandomSource = Accord.Math.Random.Generator.Random;
        }

        public void UpdateProposalDistribution(IProbabilityDistribution<double[]> proposal)
        {
            this.proposal = proposal;
        }

        /// <summary>
        ///   Attempts to generate a new observation from the target
        ///   distribution, storing its value in the <see cref="CurrentSample"/>
        ///   property.
        /// </summary>
        /// 
        /// <param name="sample">A new observation, if the method has succeed; otherwise, null.</param>
        /// 
        /// <returns>True if the sample was successfully generated; otherwise, returns false.</returns>
        /// 
        public bool TryGenerate(out double[] sample)
        {
            if (TryGenerate())
            {
                sample = current;
                return true;
            }
            else
            {
                sample = null;
                return false;
            }
        }

        /// <summary>
        ///   Attempts to generate a new observation from the target
        ///   distribution, storing its value in the <see cref="CurrentSample"/>
        ///   property.
        /// </summary>
        /// 
        /// <returns>True if the sample was successfully generated; otherwise, false.</returns>
        /// 
        public bool TryGenerate()
        {
            candidate = CandidateSample;
            Random source = RandomSource;
            double pNext = new double();
            if (targetPDF == null)
            {
                pNext = priorPDF(candidate) * likelihoodPDF(ModelOutput);
            }
            else
            {
                candidate = proposal.Generate();
                pNext = targetPDF(candidate);
            }

            if (islogPDF == false)
            {
                double Ratio = pNext / pCurrent;
                steps++;

                if (source.NextDouble() < Ratio)
                {
                    var aux = current;
                    current = candidate;
                    candidate = aux;
                    pCurrent = pNext;
                    accepts++;
                    return true;
                }

                return false;
            }
            else
            {
                double Ratio = pNext - pCurrent;
                steps++;

                if (Math.Log(source.NextDouble()) < Ratio)
                {
                    var aux = current;
                    current = candidate;
                    candidate = aux;
                    pCurrent = pNext;
                    accepts++;
                    return true;
                }

                return false;
            }
        }

        /// <summary>
        ///   Thermalizes the sample generation process, generating up to
        ///   <see cref="Discard"/> samples and discarding them. This step
        ///   is done automatically upon the first call to any of the 
        ///   <see cref="Generate()"/> functions.
        /// </summary>
        /// 
        public void WarmUp()
        {
            for (int i = 0; i < discard; i++)
                TryGenerate();
            initialized = true;
        }


        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// 
        /// <returns>
        ///   A random vector of observations drawn from this distribution.
        /// </returns>
        /// 
        public double[][] Generate(int samples)
        {
            return Generate(samples, new double[samples][]);
        }

        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// 
        /// <returns>
        ///   A random vector of observations drawn from this distribution.
        /// </returns>
        /// 
        public double[][] Generate(int samples, double[][] result)
        {
            if (!initialized)
                WarmUp();

            for (int i = 0; i < samples; i++)
            {
                while (!TryGenerate())
                {
                }
                result[i] = new double[current.Length];
                for (int j = 0; j < current.Length; j++)
                    result[i][j] = current[j];
            }

            return result;
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>
        ///   A random observation drawn from this distribution.
        /// </returns>
        /// 
        public double[] Generate()
        {
            if (!initialized)
                WarmUp();

            while (!TryGenerate()) { }

            return current;
        }

        public static MetropolisHastings<double, Independent<NormalDistribution>> Continuous(int dimensions, Func<double[], double> target)
        {
            return new MetropolisHastings<double, Independent<NormalDistribution>>(dimensions, target,
                new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()));
        }

        public static MetropolisHastings<double, Independent<NormalDistribution>, T> Continuous<T>(int dimensions, T distribution)
            where T : IMultivariateDistribution<double[]>
        {
            return new MetropolisHastings<double, Independent<NormalDistribution>, T>(distribution,
                new Independent<NormalDistribution>(dimensions, () => new NormalDistribution()));
        }

        public static MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>> Discrete(int dimensions, Func<int[], double> target)
        {
            return new MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>>(dimensions, target,
                new Independent<SymmetricGeometricDistribution, int>(dimensions, () => new SymmetricGeometricDistribution(0.5)));
        }

        public static MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>, T> Discrete<T>(int dimensions, T distribution)
            where T : IMultivariateDistribution<int[]>
        {
            return new MetropolisHastings<int, Independent<SymmetricGeometricDistribution, int>, T>(distribution,
                new Independent<SymmetricGeometricDistribution, int>(dimensions, () => new SymmetricGeometricDistribution(0.5)));
        }

    }

}

