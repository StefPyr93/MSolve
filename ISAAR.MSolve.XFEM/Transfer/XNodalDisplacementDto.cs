using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Transfer
{
    [Serializable]
    public class XNodalDisplacementDto
    {
        public double amount;
        public int node;
        public int dofID;

        public XNodalDisplacementDto(XNode node, Constraint constraint, IDofSerializer dofSerializer)
        {
            this.node = node.ID;
            this.amount = constraint.Amount;
            this.dofID = dofSerializer.Serialize(constraint.DOF);
        }

        public void Deserialize(IReadOnlyDictionary<int, XNode> allNodes, IDofSerializer dofSerializer)
        {
            XNode targetNode = allNodes[this.node];
            var constraint = new Constraint() { Amount = this.amount, DOF = dofSerializer.Deserialize(this.dofID) };
            targetNode.Constraints.Add(constraint);
        }
    }
}
