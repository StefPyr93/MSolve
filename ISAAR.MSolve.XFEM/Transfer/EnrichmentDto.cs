//using System;
//using System.Collections.Generic;
//using System.Text;
//using ISAAR.MSolve.Discretization.Transfer;
//using ISAAR.MSolve.FEM.Interpolation;
//using ISAAR.MSolve.Geometry.Coordinates;
//using ISAAR.MSolve.XFEM.Elements;
//using ISAAR.MSolve.XFEM.Enrichments.Items;
//using ISAAR.MSolve.XFEM.Entities;
//using ISAAR.MSolve.XFEM.FreedomDegrees;
//using ISAAR.MSolve.XFEM.Utilities;

////TODO: For now the crack description will be contained only in master process. It will evaluate the enrichments and store it in
////      XNodes, which will then be transfered to other processes. This class serves to transfer only the data needed for 
////      connectivity related operations. FUCK!!! EvaluateFunctionsAt() will be needed by XContinuumElement and it needs LSM data.
//namespace ISAAR.MSolve.XFEM.Transfer
//{
//    public class EnrichmentDto 
//    {
//        public int[] dofIds;

//        public EnrichmentDto(IEnrichmentItem2D enrichment, IDofSerializer dofSerializer)
//        {
//            this.dofIds = new int[enrichment.Dofs.Count];
//            for (int i = 0; i < dofIds.Length; ++i) dofIds[i] = dofSerializer.Serialize(enrichment.Dofs[i]);
//        }

//        public IEnrichmentItem2D Deserialize(IDofSerializer dofSerializer)
//        {
//            var enrichedDofs = new EnrichedDof[dofIds.Length];
//            for (int i = 0; i < dofIds.Length; ++i) enrichedDofs[i] = dofSerializer.Deserialize(dofIds[i]);
//            return new EnrichmentItemPlaceholder(enrichedDofs);
//        }

//        private class EnrichmentItemPlaceholder : IEnrichmentItem2D
//        {
//            public EnrichmentItemPlaceholder(EnrichedDof[] dofs)
//            {
//                this.Dofs = dofs;
//            }

//            public IReadOnlyList<EnrichedDof> Dofs { get; }

//            public EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint point, XContinuumElement2D element,
//                EvalInterpolation2D interpolation)
//            {
//                throw new NotImplementedException();
//            }

//            public double[] EvaluateFunctionsAt(XNode node)
//            {
//                throw new NotSupportedException("This object is only used to transfer basic data to other processes."
//                    + " For functionality only the original one can be used.");
//            }

//            public IReadOnlyList<CartesianPoint> IntersectionPointsForIntegration(XContinuumElement2D element)
//            {
//                throw new NotSupportedException("This object is only used to transfer basic data to other processes."
//                    + " For functionality only the original one can be used.");
//            }
//        }
//    }
//}
