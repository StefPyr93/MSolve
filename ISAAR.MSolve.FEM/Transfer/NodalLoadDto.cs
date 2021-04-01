using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Transfer
{
    [Serializable]
    public class NodalLoadDto
    {
        public double amount;
        public int node;
        public int dofID;

        public NodalLoadDto(Load load, IDofSerializer dofSerializer)
        {
            this.node = load.Node.ID;
            this.dofID = dofSerializer.Serialize(load.DOF);
            this.amount = load.Amount;
        }

        public Load Deserialize(IReadOnlyDictionary<int, Node> allNodes, IDofSerializer dofSerializer)
            => new Load() { Node = allNodes[this.node], DOF = dofSerializer.Deserialize(this.dofID), Amount = this.amount };
    }
}
