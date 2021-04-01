using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.Logging;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    /// <summary>
    /// The interface problem is solved using PCG. The matrix of the coarse problem KccStar, namely the static condensation of 
    /// the remainder dofs onto the corner dofs is performed explicitly.
    /// </summary>
    public class FetiDPInterfaceProblemSolverSerial : IFetiDPInterfaceProblemSolver
    {
        private readonly IModel model;
        private readonly PcgSettings pcgSettings;

        public FetiDPInterfaceProblemSolverSerial(IModel model, PcgSettings pcgSettings)
        {
            this.model = model;
            this.pcgSettings = pcgSettings;
        }

        public Vector PreviousLambda { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public bool UsePreviousLambda { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public ReorthogonalizedPcg Pcg => throw new NotImplementedException();

        public bool UseStagnationCriterion { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public Vector SolveInterfaceProblem(IFetiDPMatrixManager matrixManager,
            ILagrangeMultipliersEnumerator lagrangesEnumerator, IFetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, double globalForcesNorm, ISolverLogger logger)
        {
            int systemOrder = flexibility.NumGlobalLagrangeMultipliers;

            // Prepare PCG matrix, preconditioner, rhs and solution
            var pcgMatrix = new FetiDPInterfaceProblemMatrixSerial(matrixManager, flexibility);
            var pcgPreconditioner = new FetiDPInterfaceProblemPreconditioner(preconditioner);
            Vector globalDr = CalcGlobalDr(matrixManager, lagrangesEnumerator);
            Vector pcgRhs = CalcInterfaceProblemRhs(matrixManager, flexibility, globalDr);
            var lagranges = Vector.CreateZero(systemOrder);


            


            // Solve the interface problem using PCG algorithm
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.MaxIterationsProvider = pcgSettings.MaxIterationsProvider;
            pcgBuilder.ResidualTolerance = pcgSettings.ConvergenceTolerance;
            pcgBuilder.Convergence = pcgSettings.ConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            PcgAlgorithm pcg = pcgBuilder.Build(); //TODO: perhaps use the pcg from the previous analysis if it has reorthogonalization.
            IterativeStatistics stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagranges, true,
                () => Vector.CreateZero(systemOrder));

            // Log statistics about PCG execution
            FetiDPInterfaceProblemUtilities.CheckConvergence(stats);
            logger.LogIterativeAlgorithm(stats.NumIterationsRequired, stats.ResidualNormRatioEstimation);

            #region debug1
            var cnstValues = new CnstValues();
            if (CnstValues.printInterfaceSolutionStats)
            {
                int nL = lagranges.Length;
                int nC = matrixManager.CoarseProblemRhs.Length;
                var writer = new LinearAlgebra.Output.FullMatrixWriter();


                //// Process FIrr
                Matrix FIrr = MultiplyWithIdentity(nL, nL, flexibility.MultiplyFIrr);
                FIrr = 0.5 * (FIrr + FIrr.Transpose());

                if (CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats) {
                    string pathFIrr = cnstValues.solverPath + @"\a_FIrr_fetiDP.txt";
                    string pathFIrc = cnstValues.solverPath + @"\a_FIrc_fetiDP.txt";
                    writer.WriteToFile(FIrr, pathFIrr);
                    Matrix FIrc = MultiplyWithIdentity(nL, nC, (x, y) => y.CopyFrom(flexibility.MultiplyFIrc(x)));
                    writer.WriteToFile(FIrr, pathFIrc);
                    string pathRhs = cnstValues.solverPath + @"\a_pcg_rhs_fetiDP.txt";
                    new LinearAlgebra.Output.FullVectorWriter().WriteToFile(pcgRhs, pathRhs); }

                //(Matrix rrefFIrr, List<int> independentColsFIrr) = FIrr.ReducedRowEchelonForm();

                bool isFIrrInvertible = false;
                double detFIrr = double.NaN;
                try
                {
                    detFIrr = FIrr.CalcDeterminant();
                    isFIrrInvertible = true;
                }
                catch (Exception) { }


                bool isFIrrPosDef = false;
                try
                {
                    SkylineMatrix skyFIrr = SkylineMatrix.CreateFromMatrix(FIrr);
                    double tol = 1E-50;
                    var FIrrFactorized = skyFIrr.FactorCholesky(false, tol);
                    isFIrrPosDef = true;
                }
                catch (Exception) { }


                //// Process PCG matrix
                Matrix pcgMatrixExplicit = MultiplyWithIdentity(nL, nL, pcgMatrix.Multiply);
                pcgMatrixExplicit = 0.5 * (pcgMatrixExplicit + pcgMatrixExplicit.Transpose());

                if (CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats)
                {
                    string pathPcgMatrix = cnstValues.solverPath + @"\a_pcgMAT_fetiDP.txt";
                    writer.WriteToFile(pcgMatrixExplicit, pathPcgMatrix);
                }
                //(Matrix rref, List<int> independentCols) = pcgMatrixExplicit.ReducedRowEchelonForm();


                #region debug specific node data
                /*
                var nodeCoords = new double[3] { -22.5, -22.5, -22.5 };
                int nodeId = ((Model)model).NodesDictionary.Values.Where(x => ((x.X1 == nodeCoords[0]) && (x.X2 == nodeCoords[1]) && (x.X3 == nodeCoords[2]))).ToList().ElementAt(0).ID;
                int subdID = ((Model)model).NodesDictionary[nodeId].SubdomainsDictionary.ElementAt(0).Key;
                var nodeCoords2 = new double[3] { 0, -22.5, -22.5 };
                int nodeId2 = ((Model)model).NodesDictionary.Values.Where(x => ((x.X1 == nodeCoords2[0]) && (x.X2 == nodeCoords2[1]) && (x.X3 == nodeCoords2[2]))).ToList().ElementAt(0).ID;
                var cornerNodeCoords = new double[3] { 0, 0, 0 };
                int cornerNodeId = ((Model)model).NodesDictionary.Values.Where(x => ((x.X1 == cornerNodeCoords[0]) && (x.X2 == cornerNodeCoords[1]) && (x.X3 == cornerNodeCoords[2]))).ToList().ElementAt(0).ID;
                var orddering = ((FetiDPDofSeparatorSerial)((LagrangeMultipliersEnumeratorSerial)lagrangesEnumerator).dofSeparator).GetCornerDofOrdering(model.GetSubdomain(subdID));
                int corner_dof_x = orddering[model.GetNode(cornerNodeId), StructuralDof.TranslationX];

                int node2_x_lagrange_ID = 0;// = lagrangesEnumerator. 
                for (int i1 = 0; i1 < lagrangesEnumerator.LagrangeMultipliers.Count(); i1++)
                {
                    int LagrangeNodeID = lagrangesEnumerator.LagrangeMultipliers[i1].Node.ID;
                    var LagrangeDoftype = lagrangesEnumerator.LagrangeMultipliers[i1].DofType;
                    if ((LagrangeNodeID == nodeId2) && (LagrangeDoftype == StructuralDof.TranslationX))
                    {
                        node2_x_lagrange_ID = i1;
                    }
                }

                double node2_x_pcg_mat_value = pcgMatrixExplicit[node2_x_lagrange_ID, node2_x_lagrange_ID];
                double node2_x_pcg_rhs_value = pcgRhs[node2_x_lagrange_ID];

                var crossPointCoords = new double[3][] { new double[3] { -22.5, 0, 0 }, new double[3] { 0, -22.5, 0 }, new double[3] { 0, 0, -22.5 } };
                int[] crossPointIds = crossPointCoords.Select(x => ((Model)model).NodesDictionary.Values.Where(y => ((y.X1 == x[0]) && y.X2 == x[1]) && (y.X3 == x[2])).ToList().ElementAt(0).ID).ToArray();
                List<int> cross_point_lagranges_ID = new List<int>();
                for (int i1 = 0; i1 < lagrangesEnumerator.LagrangeMultipliers.Count(); i1++)
                {
                    int LagrangeNodeID = lagrangesEnumerator.LagrangeMultipliers[i1].Node.ID;
                    var LagrangeDoftype = lagrangesEnumerator.LagrangeMultipliers[i1].DofType;
                    if ((LagrangeNodeID == crossPointIds[0]) && (LagrangeDoftype == StructuralDof.TranslationX)) //TODO 2
                    {
                        cross_point_lagranges_ID.Add(i1);
                    }
                }
                double[] crossPoint_x_pcg_mat_values = cross_point_lagranges_ID.Select(x => pcgMatrixExplicit[x, x]).ToArray();
                double[] crossPoint_x_pcg_rhs_values = cross_point_lagranges_ID.Select(x => pcgRhs[x]).ToArray(); //DIAFORETIKO
                double[] crossPoint_x_pcg_Dr_values = cross_point_lagranges_ID.Select(x => globalDr[x]).ToArray();

                double[] crossPoint_x_Firr_mat_values = cross_point_lagranges_ID.Select(x => FIrr[x, x]).ToArray();
                double[] crossPoint_x_node_2_x_Firr_mat_values = cross_point_lagranges_ID.Select(x => FIrr[x, node2_x_lagrange_ID]).ToArray();
                double node_2_x_Firr_values = FIrr[node2_x_lagrange_ID, node2_x_lagrange_ID];
                double[] crossPoint_x_Firc_rhs_values = cross_point_lagranges_ID.Select(x => FIrc[x, corner_dof_x]).ToArray();

                // Use reflection to set the necessary matrices
                var ch01 = matrixManager.GetFetiDPSubdomainMatrixManager(model.GetSubdomain(subdID));
                FieldInfo fi;
                fi = typeof(FetiDPSubdomainMatrixManagerSkyline).GetField("Krr", BindingFlags.NonPublic | BindingFlags.Instance);
                SkylineMatrix Krr =(SkylineMatrix)(fi.GetValue(ch01));

                fi = typeof(FetiDPSubdomainMatrixManagerSkyline).GetField("Krc", BindingFlags.NonPublic | BindingFlags.Instance);
                CscMatrix Krc = (CscMatrix)(fi.GetValue(ch01));

                var dofSeparator = ((FetiDPFlexibilityMatrixSerial)flexibility).dofSeparator;
                DofTable subdRemainderDofs = dofSeparator.GetRemainderDofOrdering(model.GetSubdomain(subdID));
                DofTable subdCornerDofs = dofSeparator.GetCornerDofOrdering(model.GetSubdomain(subdID));

                int cornerNode_xdof_order = subdCornerDofs[model.GetNode(cornerNodeId), StructuralDof.TranslationX];
                int node2_remainder_xdof_order = subdRemainderDofs[model.GetNode(nodeId2), StructuralDof.TranslationX];

                double krc_value = Krc[node2_remainder_xdof_order, cornerNode_xdof_order];
                double krr_value = Krr[node2_remainder_xdof_order, node2_remainder_xdof_order]; */
                #endregion





                bool isPcgMatrixInvertible = false;
                double detPcgMatrix = double.NaN;
                try
                {
                    detPcgMatrix = pcgMatrixExplicit.CalcDeterminant();
                    isPcgMatrixInvertible = true;
                }
                catch (Exception) { }

                bool isPcgMatrixPosDef = false;
                try
                {
                    SkylineMatrix skyPcgMatrix = SkylineMatrix.CreateFromMatrix(pcgMatrixExplicit);
                    double tol = 1E-50;
                    var pcgMatrixFactorized = skyPcgMatrix.FactorCholesky(false, tol);
                    isPcgMatrixPosDef = true;
                }
                catch (Exception) { }


                int nIter = stats.NumIterationsRequired;

                //// Lagranges from LU
                Debug.WriteLine($"PCG total Iterations number = {nIter}");
                var lagrangesDirect = Vector.CreateZero(nL);
                pcgMatrixExplicit.FactorLU(false).SolveLinearSystem(pcgRhs, lagrangesDirect);
                double errorLagranges = (lagranges - lagrangesDirect).Norm2() / lagrangesDirect.Norm2();

                if (CnstValues.printPcgMatRhsEtc_AndInterfaceProblemStats)
                {
                    string pathErrorLagranges = cnstValues.solverPath + @"\a_errorLagranges_LU_iters_fetiDP.txt";
                    new LinearAlgebra.Output.FullVectorWriter().WriteToFile(Vector.CreateFromArray(new double[] { errorLagranges, (double)nIter }), pathErrorLagranges);
                }

                if (CnstValues.printInterfaceSolutionStats)
                {
                    PrintInterfaceSolverStats(nC, nL, nIter, errorLagranges,
                        isFIrrInvertible, isFIrrPosDef, isPcgMatrixInvertible, isPcgMatrixPosDef);
                }
            }
            //return lagrangesDirect;
            #endregion

            return lagranges;
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
            var statsOutputPath = cnstVal.interfaceSolverStatsPath + @"\interfaceSolver_FetiDP_stats.txt";
            cnstVal.WriteToFileStringArray(statsLines, statsOutputPath);

        }
        #endregion

        private Vector CalcInterfaceProblemRhs(IFetiDPMatrixManager matrixManager, IFetiDPFlexibilityMatrix flexibility,
            Vector globalDr)
        {
            // rhs = dr - FIrc * inv(KccStar) * fcStar
            Vector temp = matrixManager.MultiplyInverseCoarseProblemMatrix(matrixManager.CoarseProblemRhs);
            temp = flexibility.MultiplyFIrc(temp);
            return globalDr - temp;
        }

        private Vector CalcGlobalDr(IFetiDPMatrixManager matrixManager, ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            var globalDr = Vector.CreateZero(lagrangesEnumerator.NumLagrangeMultipliers);
            foreach (ISubdomain sub in model.EnumerateSubdomains())
            {
                Vector subdomainDr = FetiDPInterfaceProblemUtilities.CalcSubdomainDr(sub, matrixManager, lagrangesEnumerator);
                globalDr.AddIntoThis(subdomainDr);
            }
            return globalDr;
        }
    }
}
