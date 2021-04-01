using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.LevelSetUpdating;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit
{
    /// <summary>
    /// All model and level set data exist in all processes.
    /// </summary>
    public class TrackingExteriorCrackLsmMpiRedundant : ICrackDescriptionMpi, IExteriorCrack
    {
        private readonly TrackingExteriorCrackLsm embeddedCrack;
        private readonly ProcessDistribution procs;

        public TrackingExteriorCrackLsmMpiRedundant(ProcessDistribution procs, Func<TrackingExteriorCrackLsm> createCrack)
        {
            this.procs = procs;
            this.embeddedCrack = createCrack();
        }

        public IReadOnlyList<CartesianPoint> CrackPath => embeddedCrack.CrackPath;

        public CrackBodyEnrichment2D CrackBodyEnrichment => embeddedCrack.CrackBodyEnrichment;
        public CrackTipEnrichments2D CrackTipEnrichments => embeddedCrack.CrackTipEnrichments;

        public IHeavisideSingularityResolver SingularityResolver => embeddedCrack.SingularityResolver;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll 
            => embeddedCrack.CrackBodyNodesAll;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew 
            => embeddedCrack.CrackBodyNodesNew;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified
         => embeddedCrack.CrackBodyNodesModified;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified
            => embeddedCrack.CrackBodyNodesNearModified;

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected
            =>  embeddedCrack.CrackBodyNodesRejected;

        public IReadOnlyList<CartesianPoint> CrackTips => embeddedCrack.CrackTips;

        public IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements
            => embeddedCrack.CrackTipElements;

        public IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators
            =>  embeddedCrack.CrackTipPropagators;

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew
            =>  embeddedCrack.CrackTipNodesNew;

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld
            => embeddedCrack.CrackTipNodesOld;

        public ISet<XContinuumElement2D> ElementsModified => embeddedCrack.ElementsModified;

        public IReadOnlyList<IEnrichmentItem2D> Enrichments
            => new IEnrichmentItem2D[] { CrackBodyEnrichment, CrackTipEnrichments };

        public SingleCrackLsm LevelSets { get; private set; }

        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh => embeddedCrack.Mesh;

        public IReadOnlyList<ISingleCrack> SingleCracks => new ISingleCrack[] { this };

        public SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element)
            => embeddedCrack.FindTriangleVertices(element);

        public void InitializeGeometry(CartesianPoint crackMouth, CartesianPoint crackTip)
            => embeddedCrack.InitializeGeometry(crackMouth, crackTip);

        public void InitializeGeometry(PolyLine2D initialCrack)
            => embeddedCrack.InitializeGeometry(initialCrack);

        public void Propagate(Dictionary<int, Vector> totalFreeDisplacements)
        {
            // Calculate the propagation angle and length in master process only
            double growthAngle = double.NaN;
            double growthLength = double.NaN;
            if (procs.IsMasterProcess)
            {
                (growthAngle, growthLength) = embeddedCrack.FindPropagationSegment(totalFreeDisplacements);
            }

            // Broadcast them
            procs.Communicator.Broadcast(ref growthAngle, procs.MasterProcess);
            procs.Communicator.Broadcast(ref growthLength, procs.MasterProcess);

            // Update the crack in each process
            embeddedCrack.UpdateGeometry(growthAngle, growthLength);
        }

        public void ScatterCrackData(IXModelMpi model)
        {
            // Do nothing
        }

        public double SignedDistanceOf(XNode node) => LevelSets.SignedDistanceOf(node);

        public double SignedDistanceOf(NaturalPoint point, XContinuumElement2D element, EvalInterpolation2D interpolation)
            => LevelSets.SignedDistanceOf(point, element, interpolation);

        public void UpdateEnrichments()
            => embeddedCrack.UpdateEnrichments();

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
            => embeddedCrack.UpdateGeometry(localGrowthAngle, growthLength);
    }
}
