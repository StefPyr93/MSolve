using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;

//TODO: This class should not exist. Instead the corresponding one in FETI-DP 2D should be used. The only difference is in the 
//      calculation of dr, fcStarTilde, which are assigned to MatrixManager component. For now only the 3D MatrixManager provides
//      these vectors, but they should be in 2D as well. This would also remove the need for casting
//TODO: IAugmentationConstraints augmentationConstraints are injected into the constructor since they do not exist in 2D FETI-DP.
//      Perhaps LagrangeEnumerator and MatrixManager should also be injected
//TODO: Reorder both corner and augmented dofs and store the coarse problem dof ordering. Do this in matrix manager. Then use it 
//      here (at least) for fcStarTilde. The current approach is to list all augmented dofs after all corner and apply permutations
//      before and after solving linear systems with the coarse problem matrix KccStarTilde.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// The interface problem is solved using PCG. The matrix of the coarse problem KccStar, namely the static condensation of 
    /// the remainder dofs onto the corner dofs is performed explicitly.
    /// </summary>
    public class FetiDP3dInterfaceProblemSolverSerial : IFetiDPInterfaceProblemSolver
    {
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IModel model;
        private readonly PcgSettings pcgSettings;

        public FetiDP3dInterfaceProblemSolverSerial(IModel model, PcgSettings pcgSettings,
            IAugmentationConstraints augmentationConstraints)
        {
            this.model = model;
            this.pcgSettings = pcgSettings;
            this.augmentationConstraints = augmentationConstraints;
        }

        public Vector PreviousLambda { get; set; }

        public bool UsePreviousLambda { get; set; }

        public ReorthogonalizedPcg Pcg => null;

        public bool UseStagnationCriterion { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public Vector SolveInterfaceProblem(IFetiDPMatrixManager matrixManager,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPFlexibilityMatrix flexibility,
            IFetiPreconditioner preconditioner, double globalForcesNorm, ISolverLogger logger)
        {
            int systemOrder = flexibility.NumGlobalLagrangeMultipliers;

            // Prepare PCG matrix, preconditioner, rhs and solution
            var pcgMatrix = new FetiDPInterfaceProblemMatrixSerial(matrixManager, flexibility);
            var pcgPreconditioner = new FetiDPInterfaceProblemPreconditioner(preconditioner);
            Vector globalDr = ((IFetiDP3dMatrixManager)matrixManager).GlobalDr;
            Vector pcgRhs = CalcInterfaceProblemRhs(matrixManager, flexibility, globalDr);

            Vector lagranges;
            if (!(PreviousLambda == null))
            {
                lagranges = PreviousLambda;
            }
            else
            {
                lagranges = Vector.CreateZero(systemOrder);
            }

            

            // Solve the interface problem using PCG algorithm
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.MaxIterationsProvider = pcgSettings.MaxIterationsProvider;
            pcgBuilder.ResidualTolerance = pcgSettings.ConvergenceTolerance;
            pcgBuilder.Convergence = pcgSettings.ConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            PcgAlgorithm pcg = pcgBuilder.Build(); //TODO: perhaps use the pcg from the previous analysis if it has reorthogonalization.

            IterativeStatistics stats;
            if (!(PreviousLambda == null))
            {
                stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagranges, false,
                  () => Vector.CreateZero(systemOrder));
            }
            else
            {
                stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagranges, true,
                  () => Vector.CreateZero(systemOrder));
            }

            // Log statistics about PCG execution
            FetiDPInterfaceProblemUtilities.CheckConvergence(stats);
            logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);


            #region debug1
            int nL = lagranges.Length;
            int nC = matrixManager.CoarseProblemRhs.Length;
            var writer = new LinearAlgebra.Output.FullMatrixWriter();

            var cnst = new CnstValues();
            if (CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats)
            {
                string pathRhs = (new CnstValues()).solverPath + @"\a_pcg_rhs_fetiDP3D.txt";
                new LinearAlgebra.Output.FullVectorWriter().WriteToFile(pcgRhs, pathRhs);
            }
            ////LinearAlgebra.LibrarySettings.LinearAlgebraProviders = LinearAlgebra.LinearAlgebraProviderChoice.MKL;

            //// Process FIrr
            bool isFIrrInvertible = false;
            bool isFIrrPosDef = false;
            bool isPcgMatrixInvertible = false;
            bool isPcgMatrixPosDef = false;
            if (CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats)
            {
                Matrix FIrr = MultiplyWithIdentity(nL, nL, flexibility.MultiplyFIrr);
                FIrr = 0.5 * (FIrr + FIrr.Transpose());
                SkylineMatrix skyFIrr = SkylineMatrix.CreateFromMatrix(FIrr);
                string pathFIrr = (new CnstValues()).solverPath + @"\a_FIrr_fetiDP3D.txt";
                string pathFIrc = (new CnstValues()).solverPath + @"\a_FIrc_fetiDP3D.txt";
                writer.WriteToFile(FIrr, pathFIrr);
                Matrix FIrc = MultiplyWithIdentity(nL, nC, (x, y) => y.CopyFrom(flexibility.MultiplyFIrc(x)));
                writer.WriteToFile(FIrc, pathFIrc);

                double detFIrr = double.NaN;
                try
                {
                    detFIrr = FIrr.CalcDeterminant();
                    isFIrrInvertible = true;
                }
                catch (Exception) { }



                try
                {
                    double tol = 1E-50;
                    var FIrrFactorized = skyFIrr.FactorCholesky(false, tol);
                    isFIrrPosDef = true;
                }
                catch (Exception) { }

                //// Process PCG matrix
                Matrix pcgMatrixExplicit = MultiplyWithIdentity(nL, nL, pcgMatrix.Multiply);
                pcgMatrixExplicit = 0.5 * (pcgMatrixExplicit + pcgMatrixExplicit.Transpose());
                SkylineMatrix skyPcgMatrix = SkylineMatrix.CreateFromMatrix(pcgMatrixExplicit);
                string pathPcgMatrix = (new CnstValues()).solverPath + @"\a_pcgMAT_fetiDP3D.txt";
                writer.WriteToFile(pcgMatrixExplicit, pathPcgMatrix);
                //(Matrix rref, List<int> independentCols) = pcgMatrixExplicit.ReducedRowEchelonForm();


                double detPcgMatrix = double.NaN;
                try
                {
                    detPcgMatrix = pcgMatrixExplicit.CalcDeterminant();
                    isPcgMatrixInvertible = true;
                }
                catch (Exception) { }


                try
                {
                    double tol = 1E-50;
                    var pcgMatrixFactorized = skyPcgMatrix.FactorCholesky(false, tol);
                    isPcgMatrixPosDef = true;
                }
                catch (Exception) { }

                int nIter = stats.NumIterationsRequired;

                // Lagranges from LU
                var lagrangesDirect = Vector.CreateZero(nL);
                pcgMatrixExplicit.FactorLU(false).SolveLinearSystem(pcgRhs, lagrangesDirect);
                double errorLagranges = (lagranges - lagrangesDirect).Norm2() / lagrangesDirect.Norm2();
                string pathErrorLagranges = (new CnstValues()).solverPath + @"\a_errorLagranges_LU_iters_fetiDP3D.txt";
                new LinearAlgebra.Output.FullVectorWriter().WriteToFile(Vector.CreateFromArray(new double[] { errorLagranges, (double)nIter }), pathErrorLagranges);

                if (CnstValues.printInterfaceSolutionStats)
                {
                    PrintInterfaceSolverStats(nC, nL, nIter, errorLagranges,
    isFIrrInvertible, isFIrrPosDef, isPcgMatrixInvertible, isPcgMatrixPosDef);
                }

                //Vector resDirect = pcgRhs - pcgMatrixExplicit * lagrangesDirect;
                //double normResDirect = resDirect.Norm2();

                //Vector resPcg = pcgRhs - pcgMatrixExplicit * lagranges;
                //double normResPcg = resPcg.Norm2();

                //return lagrangesDirect;
            }

            #endregion

            if (CnstValues.printInterfaceSolutionStats&&(!CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats))
            {
                PrintInterfaceSolverStats(nC, nL, stats.NumIterationsRequired, 0.00001,
                    isFIrrInvertible, isFIrrPosDef, isPcgMatrixInvertible, isPcgMatrixPosDef);
            }
            if (CnstValues.runOnlyHexaModel && (!CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats))
            {
                PrintInterfaceSolverStats(nC, nL, stats.NumIterationsRequired, 0.00001,
                    isFIrrInvertible, isFIrrPosDef, isPcgMatrixInvertible, isPcgMatrixPosDef);
            }


            return lagranges;
        }

        private void PrintInterfaceSolverStats(int nC, int nL, int nIter, double errorLagranges,
            bool isFIrrInvertible, bool isFIrrPosDef, bool isPcgMatrixInvertible, bool isPcgMatrixPosDef)
        {
            string[] statsLines = new string[] { "nCornerDofs=" + nC.ToString() + ",",
            "nLagranges=" + nL.ToString() + ",",
            "pcgIterations=" + nIter.ToString() + ",",
            "errorLagranges=" + errorLagranges.ToString() + ",",
            "isFIrrInvertible=" + isFIrrInvertible.ToString() + ",",
            "isFIrrPosDef=" + isFIrrInvertible.ToString() + ",",
            "isPcgMatrixInvertible=" + isPcgMatrixInvertible.ToString() + ",",
            "isPcgMatrixPosDef=" + isPcgMatrixPosDef.ToString() + ",",

            };

            var cnstVal = new CnstValues();
            string statsOutputPath;
            if(CnstValues.runOnlyHexaModel)
            { statsOutputPath = cnstVal.interfaceSolverStatsPath + @"\interfaceSolver_FetiDP_3d_only_hexa_stats.txt"; }
            else
            { statsOutputPath = cnstVal.interfaceSolverStatsPath + @"\interfaceSolver_FetiDP_3d_stats.txt"; }
            cnstVal.WriteToFileStringArray(statsLines, statsOutputPath);

        }

        private Vector CalcInterfaceProblemRhs(IFetiDPMatrixManager matrixManager, IFetiDPFlexibilityMatrix flexibility,
            Vector globalDr)
        {
            // rhs = dr - FIrcTilde * inv(KccStarTilde) * fcStarTilde
            Vector fcStarTilde = matrixManager.CoarseProblemRhs;
            Vector temp = matrixManager.MultiplyInverseCoarseProblemMatrix(fcStarTilde);
            temp = flexibility.MultiplyFIrc(temp);
            return globalDr - temp;
        }

        #region debug1 
        public static Matrix MultiplyWithIdentity(int numRows, int numCols, Action<Vector, Vector> matrixVectorMultiplication)
        {
            var result = Matrix.CreateZero(numRows, numCols);
            for (int j = 0; j < numCols; ++j)
            {
                var lhs = Vector.CreateZero(numCols);
                lhs[j] = 1.0;
                var rhs = Vector.CreateZero(numRows);
                matrixVectorMultiplication(lhs, rhs);
                result.SetSubcolumn(j, rhs);
            }
            return result;
        }
        #endregion
    }
}