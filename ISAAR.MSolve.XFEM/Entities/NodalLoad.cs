using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;

//TODO: Use the same class as FEM for this
namespace ISAAR.MSolve.XFEM.Entities
{
    public class NodalLoad : INodalLoad
    {
        INode INodalLoad.Node => Node;
        public XNode Node { get; }

        IDofType INodalLoad.DOF => DofType;
        public IDofType DofType { get; }

        public double Amount { get; }

        public NodalLoad(XNode node, IDofType dofType, double value)
        {
            this.Node = node;
            this.DofType = dofType;
            this.Amount = value;
        }
    }
}
