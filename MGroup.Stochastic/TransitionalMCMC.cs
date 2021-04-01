using System;
using System.Collections.Generic;
using System.Text;
using MGroup.Stochastic.Interfaces;
using Accord.Statistics;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using Accord.Math.Random;
using Accord.Math.Decompositions;
using Accord.Math;

namespace MGroup.Stochastic
{
    class TransitionalMCMC
    {
        int dimensions;
        int numSamples;
        IProbabilityDistribution<double[]> prior;
        IProbabilityDistribution<double[]> likelihood;
        Model model;
        double[][] samples;
        double[] samplesLikelihood;
        double[] weights;

        double p_down;
        double p_up;
        double p_current;
        double b;
        double S;

        public TransitionalMCMC(int dimensions, int numSamples, IProbabilityDistribution<double[]> prior, IProbabilityDistribution<double[]> likelihood)
        {
            this.dimensions = dimensions;
            this.numSamples = numSamples;
            this.prior = prior;
            this.likelihood = likelihood;
            Initialize();
        }

        public TransitionalMCMC(int dimensions, int numSamples, IProbabilityDistribution<double[]> prior, IProbabilityDistribution<double[]> likelihood, Model model)
        {
            this.dimensions = dimensions;
            this.numSamples = numSamples;
            this.prior = prior;
            this.likelihood = likelihood;
            this.model = model;
            Initialize();
        }

        void Initialize()
        {
            samples = new double[numSamples][];
            for (int i = 0; i < numSamples; i++)
            {
                samples[i] = new double[dimensions];
                samples[i] = prior.Generate();
                samplesLikelihood[i] = CalculateLikelihood(samples[i]);
            }
            weights = new double[numSamples];
            p_down = 0;
            p_up = 2;
            p_current = 0;
            b = 0.2;
            S = 1;
        }

        public double[][] GeneratePosteriorSamples()
        {
            while (p_current < 1)
            {
                var p_temp = new double();
                var sumWeights = new double();
                while (p_up - p_down > 0.0001)
                {
                    p_temp = (p_up + p_down) / 2;
                    for (int i = 0; i < numSamples; i++)
                    {
                        weights[i] = Math.Pow(samplesLikelihood[i], p_temp - p_current);
                        sumWeights += weights[i];
                    }
                    var CoV = Measures.StandardDeviation(weights) / Measures.Mean(weights);
                    if (CoV < 1) p_up = p_temp;
                    else p_down = p_temp;
                }
                p_current = p_temp;
                var meanWeights = sumWeights / numSamples;
                S = S * meanWeights;
                if (p_current > 1) p_current = 1;
                var meanSamples = new double[dimensions];
                var stdSamples = new double[dimensions, dimensions];
                for (int i = 0; i < numSamples; i++)
                {
                    weights[i] = weights[i] / sumWeights;
                    for (int j = 0; j < dimensions; j++)
                    meanSamples[j] =+ weights[i] * samples[i][j];
                }
                for (int i = 0; i < numSamples; i++)
                    for (int ii = 0; ii < dimensions; ii++)
                        for (int jj = 0; jj < dimensions; jj++)
                            stdSamples[ii, jj] += b * b * weights[i] * (samples[i][ii] - meanSamples[ii]) * (samples[i][jj] - meanSamples[jj]);
                var SVD = new SingularValueDecomposition(stdSamples);
                var diagMatrix = SVD.DiagonalMatrix;
                for (int i = 0; i < dimensions; i++)
                    diagMatrix[i, i] = Math.Sqrt(diagMatrix[i, i]);
                stdSamples = SVD.LeftSingularVectors.Dot(diagMatrix).DotWithTransposed(SVD.RightSingularVectors);
                int[] randomSampleInd = GeneralDiscreteDistribution.Random(weights, numSamples);
                for (int i = 0; i < numSamples; i++)
                {
                    double[] deviation = NormalDistribution.Random(0, 1, dimensions);
                    double[] candSample = new double[dimensions];
                    for (int ii = 0; ii < dimensions; ii++)
                        for (int jj = 0; jj < dimensions; jj++)
                            deviation[ii] = stdSamples[ii, jj] * deviation[ii];
                    for (int ii = 0; ii < dimensions; ii++)
                        candSample[ii] = samples[randomSampleInd[i]][ii] + deviation[ii];
                    var candLikelihood = CalculateLikelihood(candSample);
                    var likelihoodRatio = candLikelihood / samplesLikelihood[i];
                    if (likelihoodRatio > Generator.Random.NextDouble())
                    {
                        samples[i] = candSample;
                        samplesLikelihood[i] = candLikelihood;
                    }
                    else
                    {
                        samples[i] = samples[randomSampleInd[i]];
                        samplesLikelihood[i] = samplesLikelihood[randomSampleInd[i]];
                    }
                }
            }
            return samples;
        }

        double CalculateLikelihood(double[] sample)
        {
            var calculatedLikelihood = new double();
            if (model == null)
                calculatedLikelihood = likelihood.ProbabilityFunction(sample);
            else
            {
                var modelOutput = SolveModel(sample);
                calculatedLikelihood = likelihood.ProbabilityFunction(modelOutput);
            }
            return calculatedLikelihood;
        }

         double[] SolveModel(double[] sample)
        {
            var modelOutput = new double[0];
            return modelOutput;
        }
    }
}
