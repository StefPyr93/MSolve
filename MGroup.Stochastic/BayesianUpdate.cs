using System;
using System.Collections.Generic;
using System.Text;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Distributions.Multivariate;

namespace MGroup.Stochastic
{
    public class BayesianUpdate
    {
        //public Func<double[], double> PriorDistribution { get; private set; }
        //public Func<double[], double> LikelihoodFunction { get; private set; }
        public MultivariateContinuousDistribution PriorDistribution { get; private set; }
        public MultivariateContinuousDistribution LikelihoodFunction { get; private set; }
        public MultivariateNormalDistribution ProposalDistribution { get; private set; }

        double[] measurementValues;
        double measurementError;

        public BayesianUpdate(MultivariateContinuousDistribution priorDistribution, double[] measurementValues, double measurementError, MultivariateNormalDistribution proposalDistribution)
        {
            this.PriorDistribution = priorDistribution;
            this.ProposalDistribution = proposalDistribution;
            this.measurementValues = measurementValues;
            this.measurementError = measurementError;
            LikelihoodFunction = CreateLikelihoodFunction();
        }

        public BayesianUpdate(MultivariateContinuousDistribution priorDistribution, double[] measurementValues, double measurementError)
        {
            this.PriorDistribution = priorDistribution;
            this.ProposalDistribution = new MultivariateNormalDistribution(new double[PriorDistribution.Dimension]);
            this.measurementValues = measurementValues;
            this.measurementError = measurementError;
            LikelihoodFunction = CreateLikelihoodFunction();
        }

        private MultivariateContinuousDistribution CreateLikelihoodFunction()
        {
            var measurementErrors = new double[measurementValues.Length, measurementValues.Length];
            for (int i = 0; i < measurementErrors.GetLength(0); i++)
            {
                measurementErrors[i, i] = Math.Pow(measurementError,2);
            }

            return new MultivariateNormalDistribution(measurementValues, measurementErrors);
        }

        public double[][] GenerateSamples(int numSamples)
        {
            var samples = new double[numSamples][];
            Func<double[], double> priorPDF = PriorDistribution.ProbabilityDensityFunction;
            Func<double[], double> likelihoodPDF = LikelihoodFunction.ProbabilityDensityFunction;
            var MH = new MetropolisHastings(PriorDistribution.Dimension, priorPDF, likelihoodPDF);
            int i = 0;
            while (i < numSamples)
            {
                MH.CandidateSample = ProposalDistribution.Generate();
                MH.ModelOutput = new double[0];  //SolveModel(candidateSample)
                var accepted = MH.TryGenerate();
                if (accepted == true)
                {
                    samples[i] = MH.CandidateSample;
                    ProposalDistribution = new MultivariateNormalDistribution(MH.CurrentSample, ProposalDistribution.Covariance);
                    MH.UpdateProposalDistribution(ProposalDistribution);
                    i++;
                }
            }
            return samples;
        }

    }
}
