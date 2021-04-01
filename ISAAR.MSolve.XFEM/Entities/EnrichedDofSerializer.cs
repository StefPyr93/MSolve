using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Enrichments.Items;

namespace ISAAR.MSolve.XFEM.Entities
{
    public class EnrichedDofSerializer : IDofSerializer
    {
        private readonly Dictionary<IDofType, int> dofs2IDs;
        private readonly IDofType[] ids2Dofs;

        public EnrichedDofSerializer(ICrackDescription crack)
        {
            var enrichedDofs = new List<IDofType>(); //TODO: These must be distinct. Perhaps I should enforce that. 
            foreach (IEnrichmentItem2D enrichment in crack.Enrichments)
            {
                foreach (IDofType dof in enrichment.Dofs) enrichedDofs.Add(dof);
            }

            var ids2DofsList = new List<IDofType>(new IDofType[]
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ,
                StructuralDof.RotationX, StructuralDof.RotationY, StructuralDof.RotationZ,
                ThermalDof.Temperature, PorousMediaDof.Pressure
            });
            for (int i = 0; i < enrichedDofs.Count; ++i) ids2DofsList.Add(enrichedDofs[i]);
            ids2Dofs = ids2DofsList.ToArray();

            dofs2IDs = new Dictionary<IDofType, int>();
            dofs2IDs[StructuralDof.TranslationX] = 0;
            dofs2IDs[StructuralDof.TranslationY] = 1;
            dofs2IDs[StructuralDof.TranslationZ] = 2;
            dofs2IDs[StructuralDof.RotationX] = 3;
            dofs2IDs[StructuralDof.RotationY] = 4;
            dofs2IDs[StructuralDof.RotationZ] = 5;
            dofs2IDs[ThermalDof.Temperature] = 6;
            dofs2IDs[PorousMediaDof.Pressure] = 7;
            int numStdDofs = dofs2IDs.Count;
            for (int i = 0; i < enrichedDofs.Count; ++i) dofs2IDs[enrichedDofs[i]] = numStdDofs + i;
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
