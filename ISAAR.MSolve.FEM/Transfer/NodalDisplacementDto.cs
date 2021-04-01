using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Transfer
{
    [Serializable]
    public class NodalDisplacementDto
    {
        public double amount;
        public int node;
        public int dofID;

        public NodalDisplacementDto(Node node, Constraint constraint, IDofSerializer dofSerializer)
        {
            this.node = node.ID;
            this.amount = constraint.Amount;
            this.dofID = dofSerializer.Serialize(constraint.DOF);
        }

        public void Deserialize(IReadOnlyDictionary<int, Node> allNodes, IDofSerializer dofSerializer)
        {
            Node targetNode = allNodes[this.node];
            var constraint = new Constraint() { Amount = this.amount, DOF = dofSerializer.Deserialize(this.dofID) };
            targetNode.Constraints.Add(constraint);
        }
    }
}
