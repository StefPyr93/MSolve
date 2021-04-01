using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: time logging must be refactored
//TODO: perhaps this belongs to Solvers.Assemblers, since the vector type depends on the solver. In that case, the 
//      elementMatrixProvider should be injected by the problem/provider.
//TODO: This class may not work for coupled problems. Then I would need an interface and this class would be the default 
//      for non-coupled problems.
namespace ISAAR.MSolve.Analyzers.Loading
{
    /// <summary>
    /// Calculates the equivalent nodal forces (at the subdomain level) due to Dirichlet boundary conditions.
    /// Authors: Maria Tavlaki
    /// </summary>
    public class DirichletEquivalentLoadsAssembler : IDirichletEquivalentLoadsAssembler
    {
        //TODO: not sure if df = K * du is the best way to calcuate df.
        private IElementMatrixProvider elementProvider; 

        public DirichletEquivalentLoadsAssembler(IElementMatrixProvider elementProvider)
        {
            this.elementProvider = elementProvider;
        }

        //TODO: Should this method take the current solution vector as parameter and use it to calculate each element's displacement?
        //      Right now the node.Constraint of each elements nodes are used for this purpose, but I do not like that.
        public void ApplyEquivalentNodalLoads(ISubdomain subdomain, IVector rhs)
        {
            try
            {
            //var times = new Dictionary<string, TimeSpan>();
            //var totalStart = DateTime.Now;
            //times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            //times.Add("element", TimeSpan.Zero);
            //times.Add("addition", TimeSpan.Zero);

            // Make sure there is at least one non zero prescribed displacement. Otherwise go directly to catch clause
            // TODO: perhaps the non zero prescribed displacements should be stored separately!
            (INode node, IDofType dof, double displacement) = subdomain.Constraints.Find(du => du != 0.0);

                foreach (IElement element in subdomain.EnumerateElements()) //TODO: why go through all the elements? Most of them will not have Dirichlet bc.
                {
                    //var elStart = DateTime.Now;
                    IMatrix elementK = elementProvider.Matrix(element);

                    //double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                    //double[] localdSolution = subdomain.CalculateElementIcrementalConstraintDisplacements(element, constraintScalingFactor);
                    double[] localdSolution =
                        CalculateElementIncrementalConstraintDisplacements(subdomain, element);

                    var elementEquivalentForces = elementK.Multiply(localdSolution);

                    elementEquivalentForces.ScaleIntoThis(-1.0);
                    subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, elementEquivalentForces, rhs);

                    //times["addition"] += DateTime.Now - elStart;
                }

                //var totalTime = DateTime.Now - totalStart;
            }
            catch (KeyNotFoundException)
            {
                // There aren't any non zero prescribed displacements, therefore we do not have to calculate the equivalent 
                // nodal loads, which is an expensive operation (all elements are accessed, their stiffness is built, etc..)
            }
        }

        private static double[] CalculateElementIncrementalConstraintDisplacements(ISubdomain subdomain, IElement element)
        {
            var elementNodalDisplacements = new double[subdomain.FreeDofOrdering.CountElementDofs(element)];
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, subdomain.Constraints);
            return elementNodalDisplacements;
        }
    }
}
