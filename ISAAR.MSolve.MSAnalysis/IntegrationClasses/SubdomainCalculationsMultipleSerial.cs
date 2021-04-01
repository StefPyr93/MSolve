using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class SubdomainCalculationsMultipleSerial
    {
        public bool UseLambdaSolutionsKff { get; set; } = false;
        public Vector[] previousLambdaSolutionsKff;
                
        #region v2 methods
        public static double[][] CombineMultipleSubdomainsIntegrationVectorsIntoTotal(Dictionary<int, double[][]> VectorsSubdomains, IScaleTransitions scaleTransitions)
        {


            double[][] totalVectors = new double[scaleTransitions.MacroscaleVariableDimension()][]; //or VectorsSubdomains.getLength(0);
            int oneSubdomainID = VectorsSubdomains.Keys.ElementAt(0);
            for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
            {
                totalVectors[i1] = new double[VectorsSubdomains[oneSubdomainID][i1].GetLength(0)];
            }


            foreach (int subdomainID in VectorsSubdomains.Keys)
            {
                for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
                {
                    for (int i2 = 0; i2 < VectorsSubdomains[subdomainID][i1].GetLength(0); i2++)
                    {
                        totalVectors[i1][i2] += VectorsSubdomains[subdomainID][i1][i2];
                    }
                }
            }

            return totalVectors;

        }

        public static Dictionary<int, double[]> CalculateFppReactionsVectorSubdomains(Model model, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, Dictionary<int, IVector> solution, Dictionary<int, IVector> dSolution,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            Dictionary<int, double[]> FppReactionVectorSubdomains = new Dictionary<int, double[]>();

            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                FppReactionVectorSubdomains.Add(subdomain.ID, SubdomainCalculations.CalculateFppReactionsVector(subdomain, elementProvider, scaleTransitions, boundaryNodes,
                solution[subdomain.ID], dSolution[subdomain.ID], initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements));
            }

            return FppReactionVectorSubdomains;

        }

        internal static double[] CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal(Dictionary<int, double[]> fppReactionVectorSubdomains)
        {

            double[] totalVector = new double[fppReactionVectorSubdomains.ElementAt(0).Value.GetLength(0)];

            foreach (int subdomainID in fppReactionVectorSubdomains.Keys)
            {
                for (int i1 = 0; i1 < totalVector.GetLength(0); i1++)
                {
                    totalVector[i1] += fppReactionVectorSubdomains[subdomainID][i1];
                }
            }

            return totalVector;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="KfpDqSubdomains"></param>
        /// <param name="model"></param>
        /// <param name="elementProvider"></param>
        /// <param name="scaleTransitions"></param>
        /// <param name="boundaryNodes"></param>
        /// <param name="solver">
        /// <paramref name="solver"/>.<see cref="ISolver.Initialize"/> must already have been called. Also the linear system matrices must already have been set.
        /// </param>
        public Dictionary<int, double[][]> CalculateKffinverseKfpDqSubdomains(Dictionary<int, double[][]> KfpDqSubdomains, Model model, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, IFetiSolver solver)
        {
            //IReadOnlyDictionary<int, ILinearSystem> linearSystems = solver.LinearSystems;

            #region Creation of solution vectors structure
            Dictionary<int, double[][]> f2_vectorsSubdomains = new Dictionary<int, double[][]>();
            foreach (int subdomainID in KfpDqSubdomains.Keys)
            {
                f2_vectorsSubdomains.Add(subdomainID, new double[KfpDqSubdomains[subdomainID].GetLength(0)][]);
            }
            #endregion

            //#region Creation of linear systems with no RHS (first RHS value can be assigned too )
            //ILinearSystem[] seclinearSystems = new ILinearSystem[linearSystems.Count];
            //int counter = 0;
            //foreach (ILinearSystem subdomain in linearSystems.Values)
            //{
            //    seclinearSystems[counter] = new SkylineLinearSystem(subdomain.ID, new double[KfpDqSubdomains[subdomain.ID][0].GetLength(0)]);
            //    seclinearSystems[counter].Matrix = subdomain.Matrix;
            //}
            //#endregion

            //#region creation of solver
            //var secSolver = new SolverSkyline(seclinearSystems[0]);
            //secSolver.Initialize();

            //#endregion

            RedistributeKfpDqInSubdomainRHSs(KfpDqSubdomains, model, scaleTransitions, solver);


            if (UseLambdaSolutionsKff)
            { if (previousLambdaSolutionsKff == null) { previousLambdaSolutionsKff = new Vector[scaleTransitions.MacroscaleVariableDimension()]; } } //TODOGer1: pithanws axrhsto to copy

            #region Consecutively(for macroscaleVariableDimension times) Set proper right hand side. Solve. Copy solution in output vector 
            //int oneSubomainID = linearSystems.First().Value.Subdomain.ID;           //seclinearSystems[0].ID;
            for (int k = 0; k < scaleTransitions.MacroscaleVariableDimension(); k++) //KfpDqSubdomains[linearSystems[0].ID].GetLength(0)=Mac
            {
                #region Set proper RHS 
                //var globalRHS = new Vector(model.TotalDOFs); //TODO: uncoomment if globalRHS is needed for solver
                foreach (ISubdomain secSubdomain in model.EnumerateSubdomains())
                {
                    ILinearSystemMpi linearSystem = solver.GetLinearSystem(secSubdomain);
                    //linearSystem.Reset();
                    //TODOGer1: mhpws xreiazetai linearSystem.Solution.Clear ?s
                    linearSystem.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.ID][k], false);
                    
                    //secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], false);
                    //secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], true); Wste sigoura na mhn peiraxthei to double[]


                    //mappings[seclinearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == secSubdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
                    //TODO: uncoomment if globalRHS is needed for solver
                }
                #endregion

                #region Solve

                solver.usePreviousLambda = this.UseLambdaSolutionsKff;
                if (!(previousLambdaSolutionsKff == null)) { solver.previousLambda = previousLambdaSolutionsKff[k]; }

                #region solver parameters
                if (solver is IFetiSolver fetiSolver)
                {
                    if (fetiSolver.InterfaceProblemSolver.Pcg != null)
                    {
                        fetiSolver.InterfaceProblemSolver.Pcg.Clear();
                        //fetiSolver.InterfaceProblemSolver.Pcg.ReorthoCache.Clear();
                        if (k == 0)
                        {
                            fetiSolver.InterfaceProblemSolver.UseStagnationCriterion = true;
                            //fetiSolver.InterfaceProblemSolver.Pcg.ReorthoCache.Clear();
                        }
                        else
                        {
                            fetiSolver.InterfaceProblemSolver.UseStagnationCriterion = false;
                        }
                    }
                }
                #endregion

                solver.Solve();
                if (UseLambdaSolutionsKff)
                { previousLambdaSolutionsKff[k] = solver.previousLambda.Copy(); } //TODOGer1: pithanws axrhsto to copy
                #endregion

                #region Copy solution in output vector
                foreach (ISubdomain secSubdomain in model.EnumerateSubdomains())
                {
                    f2_vectorsSubdomains[secSubdomain.ID][k] = solver.GetLinearSystem(secSubdomain).Solution.CopyToArray();
                }
                #endregion
            }
            #endregion

            return f2_vectorsSubdomains;
        }

        private static void RedistributeKfpDqInSubdomainRHSs(Dictionary<int, double[][]> kfpDqSubdomains, Model model, IScaleTransitions scaleTransitions, ISolverMpi solver)
        {
            var redistributedkfpDqSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                #region Create KfpDq 
                redistributedkfpDqSubdomains[subdomain.ID] = new double[scaleTransitions.MacroscaleVariableDimension()][];
                for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                {
                    redistributedkfpDqSubdomains[subdomain.ID][j1] = new double[subdomain.FreeDofOrdering.NumFreeDofs]; //v2.2 subdomain.TotalDOFs]; 
                }
                #endregion
            }

            var freeNodes = model.GlobalDofOrdering.GlobalFreeDofs.GetRows();

            var loadDistributor = solver.NodalLoadDistributor;

            foreach(var freeNode in freeNodes)
            {
                var subdomains = freeNode.SubdomainsDictionary.Values;
                bool isNodeBoundaryRemainder = false;
                if (subdomains.Count > 1) isNodeBoundaryRemainder = true;

                if(isNodeBoundaryRemainder)
                {
                    #region Gather rhs contributions from all subdomains 
                    var dofTypes = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(freeNode).ToArray();
                    double[,] totalRhsValue = new double[dofTypes.Count(),scaleTransitions.MacroscaleVariableDimension()];
                    
                    foreach (ISubdomain subdomain in subdomains)
                    {
                        var dofRows = dofTypes.Select(x => subdomain.FreeDofOrdering.FreeDofs[freeNode, x]).ToArray();
                        //bool isFree = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(freeNode, elementDOFTypes[i][dofTypeRowToNumber],out int dofRow); 

                        for (int i1 = 0; i1 < dofRows.Count(); i1++)
                        {
                            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                            {
                                totalRhsValue[i1, j1] += kfpDqSubdomains[subdomain.ID][j1][dofRows[i1]];
                            }
                        }
                    }
                    #endregion

                    #region Distribute rhs
                    foreach (ISubdomain subdomain in subdomains)
                    {
                        var dofRows = dofTypes.Select(x => subdomain.FreeDofOrdering.FreeDofs[freeNode, x]).ToArray();
                        //bool isFree = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(freeNode, elementDOFTypes[i][dofTypeRowToNumber],out int dofRow); 

                        for (int i1 = 0; i1 < dofRows.Count(); i1++)
                        {
                            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                            {
                                //totalRhsValue[i1, j1] += kfpDqSubdomains[subdomain.ID][j1][i1];
                                kfpDqSubdomains[subdomain.ID][j1][dofRows[i1]] = loadDistributor.ScaleNodalLoad(subdomain, new Load() { Amount = totalRhsValue[i1, j1], DOF = dofTypes[i1], Node= (Node)freeNode });

                            }
                        }
                    }
                    #endregion

                }



            }


        }

        public static Dictionary<int, double[][]> CalculateKpfKffinverseKfpDqSubdomains(Dictionary<int, double[][]> f2_vectorsSubdomains, Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, double[][]> f3_vectorsSubdomains = new Dictionary<int, double[][]>();

            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                f3_vectorsSubdomains.Add(subdomain.ID, SubdomainCalculations.CalculateKpfKffinverseKfpDq(f2_vectorsSubdomains[subdomain.ID], subdomain, elementProvider, scaleTransitions, boundaryNodes));
            }
            return f3_vectorsSubdomains;
        }
        #endregion
    }
}
