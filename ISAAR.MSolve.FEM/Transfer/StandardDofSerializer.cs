using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Transfer;

//TODO: Transfering enriched dofs in XFEM is handled through the enrichments. Therefore this class should be refactored or at
//      least renamed.
namespace ISAAR.MSolve.FEM.Transfer
{
    /// <summary>
    /// Initialize this in each process.
    /// </summary>
    public class StandardDofSerializer : IDofSerializer
    {
        private readonly Dictionary<IDofType, int> dofs2IDs;
        private readonly IDofType[] ids2Dofs;

        public StandardDofSerializer()
        {
            dofs2IDs = new Dictionary<IDofType, int>();
            dofs2IDs[StructuralDof.TranslationX] = 0;
            dofs2IDs[StructuralDof.TranslationY] = 1;
            dofs2IDs[StructuralDof.TranslationZ] = 2;
            dofs2IDs[StructuralDof.RotationX] = 3;
            dofs2IDs[StructuralDof.RotationY] = 4;
            dofs2IDs[StructuralDof.RotationZ] = 5;
            dofs2IDs[ThermalDof.Temperature] = 6;
            dofs2IDs[PorousMediaDof.Pressure] = 7;

            ids2Dofs = new IDofType[]
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ,
                StructuralDof.RotationX, StructuralDof.RotationY, StructuralDof.RotationZ,
                ThermalDof.Temperature, PorousMediaDof.Pressure
            };
        }

        public IDofType Deserialize(int dofID)
        {
            Debug.Assert((dofID >= 0) && (dofID < ids2Dofs.Length), 
                $"Illegal dof id = {dofID}. Must belong to [0, {ids2Dofs.Length})");
            return ids2Dofs[dofID];
        }

        public int Serialize(IDofType dofType)
        {
            Debug.Assert(dofs2IDs.ContainsKey(dofType), $"Unknown dof type = {dofType}");
            return dofs2IDs[dofType];
        }
    }
}
