using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Remove FETI-DP code from her and serial implementations
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements
{
    public class FreeDofDisplacementsCalculatorMpi : IFreeDofDisplacementsCalculator
    {
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly ILagrangeMultipliersEnumerator lagrangesEnumerator;
        private readonly IFetiDPMatrixManager matrixManager;
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public FreeDofDisplacementsCalculatorMpi(ProcessDistribution processDistribution, IModel model, 
            IFetiDPDofSeparator dofSeparator, IFetiDPMatrixManager matrixManager, 
            ILagrangeMultipliersEnumerator lagrangesEnumerator)
        {
            this.procs = processDistribution;
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.matrixManager = matrixManager;
            this.lagrangesEnumerator = lagrangesEnumerator;
        }

        public void CalculateSubdomainDisplacements(Vector lagranges, IFetiDPFlexibilityMatrix flexibility)
        {
            procs.Communicator.BroadcastVector(ref lagranges, lagrangesEnumerator.NumLagrangeMultipliers, procs.MasterProcess); //TODO: Ideally calculate Br^T*lambda and scatter that
            
            //TODO: This will also broadcast lagranges. I should provide an overload of FIrc^T * vector that does not broadcast vIn.
            //      Actually lagranges are probably broadcasted in other components too. Perhaps they should be accessed by a dedicated class
            Vector uc = CalcCornerDisplacements(flexibility, lagranges); 
            
            procs.Communicator.BroadcastVector(ref uc, dofSeparator.NumGlobalCornerDofs, procs.MasterProcess);
            foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
            {
                ISubdomain subdomain = model.GetSubdomain(s);
                FreeDofDisplacementsCalculatorUtilities.CalcAndStoreFreeDisplacements(subdomain, dofSeparator, matrixManager,
                    lagrangesEnumerator, lagranges, uc);
            }
            
        }

        private Vector CalcCornerDisplacements(IFetiDPFlexibilityMatrix flexibility, Vector lagranges)
        {
            // uc = inv(KccStar) * (fcStar + FIrc^T * lagranges)
            Vector temp = flexibility.MultiplyFIrcTransposed(lagranges);
            if (procs.IsMasterProcess)
            {
                temp.AddIntoThis(matrixManager.CoarseProblemRhs);
                return matrixManager.MultiplyInverseCoarseProblemMatrix(temp);
            }
            else return null;
        }
    }
}
