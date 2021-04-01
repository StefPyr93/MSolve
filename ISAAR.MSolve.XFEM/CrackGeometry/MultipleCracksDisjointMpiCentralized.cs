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
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: replace the boilerplate code with one or more (or generic) Union() private methods
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    public class MultipleCracksDisjointMpiCentralized : ICrackDescriptionMpi
    {
        private readonly ProcessDistribution procs;
        private readonly MultipleCracksDisjoint multiCrack_master;
        private readonly List<TrackingExteriorCrackLsmMpiCentralized> singleCracks;

        public MultipleCracksDisjointMpiCentralized(ProcessDistribution procs, IReadOnlyList<TrackingExteriorCrackLsm> singleCracks)
        {
            this.procs = procs;
            this.singleCracks = new List<TrackingExteriorCrackLsmMpiCentralized>();
            foreach (TrackingExteriorCrackLsm crack in singleCracks)
            {
                this.singleCracks.Add(new TrackingExteriorCrackLsmMpiCentralized(procs, () => crack));
            }
            if (procs.IsMasterProcess) this.multiCrack_master = new MultipleCracksDisjoint(singleCracks);
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackBodyNodesAll;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackBodyNodesNew;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackBodyNodesModified;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackBodyNodesNearModified;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackBodyNodesRejected;
            }
        }

        public IReadOnlyList<CartesianPoint> CrackTips
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackTips;
            }
        }

        public IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackTipElements;
            }
        }

        public IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackTipPropagators;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackTipNodesNew;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.CrackTipNodesOld;
            }
        }

        public ISet<XContinuumElement2D> ElementsModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.ElementsModified;
            }
        }
        public IReadOnlyList<IEnrichmentItem2D> Enrichments
        {
            get
            {
                var enrichments = new List<IEnrichmentItem2D>();
                foreach (var crack in singleCracks) enrichments.AddRange(crack.Enrichments);
                return enrichments;
            }
        }


        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh
        {
            get
            {
                procs.CheckProcessIsMaster();
                return multiCrack_master.Mesh;
            }
        }

        public IReadOnlyList<ISingleCrack> SingleCracks => singleCracks;

        public void Propagate(Dictionary<int, Vector> totalFreeDisplacements)
        {
            procs.CheckProcessIsMaster();
            multiCrack_master.Propagate(totalFreeDisplacements);
        }

        public void ScatterCrackData(IXModelMpi model)
        {
            foreach (TrackingExteriorCrackLsmMpiCentralized crack in singleCracks) crack.ScatterCrackData(model);
        }

        public void UpdateEnrichments()
        {
            procs.CheckProcessIsMaster();
            multiCrack_master.UpdateEnrichments();
        }
    }
}
