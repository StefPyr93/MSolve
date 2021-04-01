using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Transfer
{
    [Serializable]
    public class XNodalLoadDto
    {
        public double amount;
        public int node;
        public int dofID;

        public XNodalLoadDto(NodalLoad load, IDofSerializer dofSerializer)
        {
            this.node = load.Node.ID;
            this.dofID = dofSerializer.Serialize(load.DofType);
            this.amount = load.Amount;
        }

        public NodalLoad Deserialize(IReadOnlyDictionary<int, XNode> allNodes, IDofSerializer dofSerializer)
            => new NodalLoad(allNodes[this.node], dofSerializer.Deserialize(this.dofID), this.amount);
    }
}
