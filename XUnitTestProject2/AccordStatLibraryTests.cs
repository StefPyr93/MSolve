using System;
using System.Diagnostics;
using Xunit;
using Accord.Math.Random;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using MGroup.Stochastic;
//using ISAAR.MSolve.Sto

namespace XUnitTestProject2
{
    public class AccordStatLibraryTests
    {
        [Fact]
        public void Test1()
        {
            //var mvn = new Accord.Statistics.Distributions.Multivariate.MultivariateNormalDistribution(new double[2] { 0, 0 });
            var target = new MultivariateNormalDistribution(new double[2] { 0, 0 });
            Debug.WriteLine("mytarget: Type is {0}", target.GetType());
            Func<double[], double> logtargetPDF = target.LogProbabilityDensityFunction;
            Func<double[], double> targetPDF = target.ProbabilityDensityFunction;
            var dim = target.Dimension;
            //var MHlog = new MetropolisHastings(dim, logtargetPDF); MHlog.islogPDF = true;
            var MH = new MetropolisHastings(dim, targetPDF);
            //var singleSample = MH.Generate();
            var moreSamples = MH.Generate(100);
            var acceptance = MH.AcceptanceRate;
        }

        [Fact]
        public void Test2()
        {
            var prior = new MultivariateNormalDistribution(new double[2] { 10, 10 });
            var likelihood = new MultivariateNormalDistribution(new double[2] { 20, 20 });
            var proposal = new MultivariateNormalDistribution(new double[2] { 0, 0 });
            var dim = prior.Dimension;
            var numSamples = 10000;
            var samples = new double[numSamples][];
            Func<double[], double> priorPDF = prior.ProbabilityDensityFunction;
            Func<double[], double> likelihoodPDF = likelihood.ProbabilityDensityFunction;
            var MH = new MetropolisHastings(dim, priorPDF, likelihoodPDF);
            int i = 0;
            while (i < numSamples)
            {
                MH.CandidateSample = proposal.Generate();
                MH.ModelOutput = new double[2] { 0, 0 };  //SolveModel(candidateSample)
                var accepted = MH.TryGenerate();
                if (accepted == true)
                {
                    samples[i] = MH.CandidateSample;
                    proposal = new MultivariateNormalDistribution(MH.CurrentSample, proposal.Covariance);
                    MH.UpdateProposalDistribution(proposal);
                    i++;
                }
            }
            var posterior = new MultivariateNormalDistribution(dim);
            posterior.Fit(samples);
        }

        [Fact]
        public void Test3()
        {
            //var probabilities = new double[3] { Math.Log(0.95), Math.Log(0.025), Math.Log(0.025) };
            var probabilities = new double[3] { 0.025, 0.95, 0.025 };
            var numberOfTrials = 5;
            var discrDistr = new GeneralDiscreteDistribution(probabilities);
            //double[] mean = discrDistr.Mean;     // {  1.25, 3.75 }
            //double[] median = discrDistr.Median; // {  1.25, 3.75 }
            //double[] var = discrDistr.Variance;  // { -0.9375, -0.9375 }
            var source = Generator.Random;
            var sample = discrDistr.Generate(source);
            var random = GeneralDiscreteDistribution.Random(probabilities,10);
        }
    }
}
