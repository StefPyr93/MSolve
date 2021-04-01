//using System;
//using System.Collections.Generic;
//using System.Text;
//using ISAAR.MSolve.XFEM.Enrichments.Items;

////TODO: Better to have a distributed LSM that will recalculate nodal level sets and enrichments
//namespace ISAAR.MSolve.XFEM.Transfer
//{
//    public class EnrichmentSerializer
//    {
//        private int numEnrichments = 0;
//        private Dictionary<IEnrichmentItem2D, int> enrichmentsToIds = new Dictionary<IEnrichmentItem2D, int>();
//        private Dictionary<int, IEnrichmentItem2D> idsToEnrichments = new Dictionary<int, IEnrichmentItem2D>();

//        public EnrichmentSerializer()
//        {
//        }

//        public void AddEnrichment(IEnrichmentItem2D enrichment)
//        {
//            enrichmentsToIds[enrichment] = numEnrichments;
//            idsToEnrichments[numEnrichments] = enrichment;
//            ++numEnrichments;
//        }

//        public IEnrichmentItem2D GetEnrichment(int id) => idsToEnrichments[id];
//        public int GetEnrichmentID(IEnrichmentItem2D enrichment) => enrichmentsToIds[enrichment];
//    }
//}
