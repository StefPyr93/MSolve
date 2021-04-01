using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: there are algebraic expressions for these. E.g. inv(Lb^T * Db * Lb) for the inverse of total stiffness. Should I use 
//      those instead?
//TODO: Could this be used outside FETI?
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class BoundaryDofLumpedStiffness
    {
        internal BoundaryDofLumpedStiffness(Dictionary<ISubdomain, double> subdomainStiffnesses, double totalStiffness)
        {
            this.SubdomainStiffnesses = subdomainStiffnesses;
            this.TotalStiffness = totalStiffness;
        }

        internal Dictionary<ISubdomain, double> SubdomainStiffnesses { get; }
        internal double TotalStiffness { get; }

        internal double CalcRelativeStiffness(ISubdomain subdomain) => SubdomainStiffnesses[subdomain] / TotalStiffness;
    }
}
