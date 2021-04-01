using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution
{
    public class FetiDPHomogeneousStiffnessDistributionOLD : HomogeneousStiffnessDistributionOLD
    {
        private readonly FetiDPDofSeparatorOLD dofSeparator;

        public FetiDPHomogeneousStiffnessDistributionOLD(IModel model, FetiDPDofSeparatorOLD dofSeparator) : base(model, dofSeparator)
        {
            this.dofSeparator = dofSeparator;
        }

        public override double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load)
        {
            INode node = load.Node;
            IDofType dof = load.DOF;

            // Loads at corner dofs will be distributed equally. It shouldn't matter how I distribute these, since I 
            // will only sum them together again during the static condensation of remainder dofs phase.
            //TODO: is that correct?
            bool isCornerDof = dofSeparator.GlobalCornerDofOrdering.Contains(node, dof);
            if (isCornerDof) return load.Amount / node.Multiplicity;
            else if (node.Multiplicity > 1) return load.Amount / node.Multiplicity;
            else return load.Amount;
        }

        //public Dictionary<int, SparseVector> DistributeNodalLoadsOLD(Dictionary<int, ISubdomain> subdomains,
        //    Table<INode, IDofType, double> globalNodalLoads)
        //    => FetiDPStiffnessDistributionUtilitiesOLD.DistributeNodalLoadsOLD(dofSeparator, subdomains, globalNodalLoads,
        //        base.CalcBoundaryDofCoefficientsOLD);
    }
}
