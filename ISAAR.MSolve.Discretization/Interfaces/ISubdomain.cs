using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: tidy up the methods that concern material state
namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain
    {
        void CalculateStressesOnly(IVectorView solution, IVectorView dSolution);

        IVector CalculateRHSonly(IVectorView solution, IVectorView dSolution);
        Table<INode, IDofType, double> Constraints { get; }

        ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        ISubdomainFreeDofOrdering FreeDofOrdering { get; set; } //TODO: this should not be managed by the subdomain. Update after 6 months: yeap, see the mess in collocation

        Vector Forces { get; set; } //TODO: this should be a Vector or IVector and stored elsewhere.

        int ID { get; }

        bool ConnectivityModified { get; set; }
        bool StiffnessModified { get; set; }

        int NumElements { get; }
        int NumNodalLoads { get; }
        int NumNodes { get; }

        void ClearMaterialStresses();

        void ConnectDataStructures();

        IEnumerable<IElement> EnumerateElements();
        IEnumerable<INodalLoad> EnumerateNodalLoads();
        IEnumerable<INode> EnumerateNodes();

        IElement GetElement(int elementID);
        INode GetNode(int nodeID);

        IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution); //TODO: this should be done by a dedicated class instead of the subdomain

        void ResetMaterialsModifiedProperty();

        void ScaleConstraints(double scalingFactor); //TODO: this should be done by a dedicated class instead of the subdomain

        void SaveMaterialState();
    }
}
