using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

//TODO: replace the boilerplate code with one or more (or generic) Union() private methods
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    public class MultipleCracksDisjointMpiRedundant : ICrackDescriptionMpi
    {
        private readonly ProcessDistribution procs;
        private readonly MultipleCracksDisjoint multiCrack;
        private readonly List<TrackingExteriorCrackLsmMpiRedundant> singleCracks;

        public MultipleCracksDisjointMpiRedundant(ProcessDistribution procs, 
            IReadOnlyList<TrackingExteriorCrackLsm> singleCracks)
        {
            this.procs = procs;
            this.singleCracks = new List<TrackingExteriorCrackLsmMpiRedundant>();
            foreach (TrackingExteriorCrackLsm crack in singleCracks)
            {
                this.singleCracks.Add(new TrackingExteriorCrackLsmMpiRedundant(procs, () => crack));
            }
            this.multiCrack = new MultipleCracksDisjoint(singleCracks);
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll 
            => multiCrack.CrackBodyNodesAll;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew
            => multiCrack.CrackBodyNodesNew;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified
            => multiCrack.CrackBodyNodesModified;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified
            => multiCrack.CrackBodyNodesNearModified;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected
            => multiCrack.CrackBodyNodesRejected;
        
        public IReadOnlyList<CartesianPoint> CrackTips => multiCrack.CrackTips;

        public IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements
            => multiCrack.CrackTipElements;

        public IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators
            => multiCrack.CrackTipPropagators;

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew
            => multiCrack.CrackTipNodesNew;

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld
            => multiCrack.CrackTipNodesOld;

        public ISet<XContinuumElement2D> ElementsModified => multiCrack.ElementsModified;

        public IReadOnlyList<IEnrichmentItem2D> Enrichments
        {
            get
            {
                var enrichments = new List<IEnrichmentItem2D>();
                foreach (var crack in singleCracks) enrichments.AddRange(crack.Enrichments);
                return enrichments;
            }
        }


        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh => multiCrack.Mesh;

        public IReadOnlyList<ISingleCrack> SingleCracks => singleCracks;

        public void Propagate(Dictionary<int, Vector> totalFreeDisplacements)
        {
            foreach (TrackingExteriorCrackLsmMpiRedundant crack in singleCracks) crack.Propagate(totalFreeDisplacements);
        }

        public void ScatterCrackData(IXModelMpi model)
        {
            // Do nothing
        }

        public void UpdateEnrichments() => multiCrack.UpdateEnrichments();
    }
}
