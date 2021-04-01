using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    public class Feti1HeterogeneousStiffnessDistribution : HeterogeneousStiffnessDistributionOLD
    {
        public Feti1HeterogeneousStiffnessDistribution(IModel model, Feti1DofSeparator dofSeparator) : base(model, dofSeparator)
        {
        }

        public override double ScaleNodalLoad(ISubdomain subdomain, INodalLoad load)
        {
            INode node = load.Node;
            if (node.Multiplicity == 1) return load.Amount;
            else return boundaryDofStiffnesses[node, load.DOF].CalcRelativeStiffness(subdomain) * load.Amount;
        }

        //public Dictionary<int, SparseVector> DistributeNodalLoadsOLD(Dictionary<int, ISubdomain> subdomains,
        //    Table<INode, IDofType, double> globalNodalLoads)
        //    => Feti1StiffnessDistributionUtilitiesOLD.DistributeNodalLoadsOLD(dofSeparator, subdomains, globalNodalLoads,
        //        base.CalcBoundaryDofCoefficientsOLD);
    }
}
