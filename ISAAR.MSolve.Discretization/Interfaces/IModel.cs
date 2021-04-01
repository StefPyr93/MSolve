using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IModel
    {
        Table<INode, IDofType, double> Constraints { get; }

        IDofSerializer DofSerializer { get; set; }

        IGlobalFreeDofOrdering GlobalDofOrdering { get; set; } //TODO: this should not be managed by the model. Update after 6 months: yeap, see the mess in collocation
        IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; }

        int NumElements { get; }
        int NumNodes { get; }
        int NumSubdomains { get; }

        //TODO: Applying loads is not the job of the model. An analyzer should decide when that happens and the job should be 
        //      delegated to dedicated classes for each load type. Distributing the loads between subdomains is then performed 
        //      by the solver, which the model has no knowledge of, but the analyzer does. Assigning the loads to the subdomains
        //      or otherwise converting them to a usable form should be incorporated in ConnectDataStructures(). This logic has 
        //      already been implemented for nodal loads and prescribed displacements. The rest should follow.
        void ApplyLoads(); //TODOMaria: Here is where the element loads are assembled
        void ApplyMassAccelerationHistoryLoads(int timeStep);

        void ConnectDataStructures();

        IEnumerable<IElement> EnumerateElements(); //TODO: At some point I must do the same for concrete classes
        IEnumerable<INode> EnumerateNodes();
        IEnumerable<ISubdomain> EnumerateSubdomains();

        IElement GetElement(int elementID);
        INode GetNode(int nodeID);
        ISubdomain GetSubdomain(int subdomainID);

        ////TODO: This circumvents the covariance issue between Dictionary<int, Node> and Dictionary<int, INode>. Is there a more elegant solution?
        ////TODO: Similarly there should be methods and properties NumNodes, EnumerateNodes(), AddNode(), RemoveNode(). 
        /////This is much better than exposing lists or dictionaries. Ditto for elements and subdomains.
        //INode GetNode(int nodeID); 
    }
}
