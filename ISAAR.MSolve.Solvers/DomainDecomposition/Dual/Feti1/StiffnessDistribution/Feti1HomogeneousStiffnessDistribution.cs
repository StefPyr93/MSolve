using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    public class Feti1HomogeneousStiffnessDistribution : HomogeneousStiffnessDistributionOLD
    {
        public Feti1HomogeneousStiffnessDistribution(IModel model, Feti1DofSeparator dofSeparator) : base(model, dofSeparator)
        {
        }

        public override double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load)
        {
            INode node = load.Node;
            if (node.Multiplicity > 1) return load.Amount / node.Multiplicity;
            else return load.Amount;
        }

        //public Dictionary<int, SparseVector> DistributeNodalLoadsOLD(Dictionary<int, ISubdomain> subdomains,
        //    Table<INode, IDofType, double> globalNodalLoads)
        //    => Feti1StiffnessDistributionUtilitiesOLD.DistributeNodalLoadsOLD(dofSeparator, subdomains, globalNodalLoads,
        //        base.CalcBoundaryDofCoefficientsOLD);
    }
}
