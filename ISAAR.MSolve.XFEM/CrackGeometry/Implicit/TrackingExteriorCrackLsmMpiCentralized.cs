using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
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
    /// Warning: may misclassify elements as tip elements, causing gross errors.
    /// </summary>
    public class TrackingExteriorCrackLsmMpiCentralized : ICrackDescriptionMpi, IExteriorCrack
    {
        private const int levelSetDataTag = 0;
        private const int enrichmentDataTag = 1;

        private readonly TrackingExteriorCrackLsm embeddedCrack_master;
        private readonly ProcessDistribution procs;

        public TrackingExteriorCrackLsmMpiCentralized(ProcessDistribution procs, Func<TrackingExteriorCrackLsm> createCrack)
        {
            this.procs = procs;
            if (procs.IsMasterProcess)
            {
                this.embeddedCrack_master = createCrack();
                this.CrackBodyEnrichment = embeddedCrack_master.CrackBodyEnrichment;
                this.CrackTipEnrichments = embeddedCrack_master.CrackTipEnrichments;
            }
            else
            {
                this.CrackBodyEnrichment = new CrackBodyEnrichment2D(this);
                this.CrackTipEnrichments = new CrackTipEnrichments2D(this, CrackTipPosition.Single);
            }
        }

        public IReadOnlyList<CartesianPoint> CrackPath
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackPath;
            }
        }

        public CrackBodyEnrichment2D CrackBodyEnrichment { get; }
        public CrackTipEnrichments2D CrackTipEnrichments { get; }

        public IHeavisideSingularityResolver SingularityResolver
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.SingularityResolver;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackBodyNodesAll;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackBodyNodesNew;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackBodyNodesModified;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackBodyNodesNearModified;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackBodyNodesRejected;
            }
        }

        public IReadOnlyList<CartesianPoint> CrackTips
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackTips;
            }
        }

        public IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackTipElements;
            }
        }

        public IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackTipPropagators;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackTipNodesNew;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.CrackTipNodesOld;
            }
        }

        public ISet<XContinuumElement2D> ElementsModified
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.ElementsModified;
            }
        }

        public IReadOnlyList<IEnrichmentItem2D> Enrichments
            => new IEnrichmentItem2D[] { CrackBodyEnrichment, CrackTipEnrichments };

        public SingleCrackLsm LevelSets { get; private set; }

        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh
        {
            get
            {
                procs.CheckProcessIsMaster();
                return embeddedCrack_master.Mesh;
            }
        }

        public IReadOnlyList<ISingleCrack> SingleCracks => new ISingleCrack[] { this };

        public SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element)
        {
            return embeddedCrack_master.FindTriangleVertices(element);
        }

        public void InitializeGeometry(CartesianPoint crackMouth, CartesianPoint crackTip)
        {
            procs.CheckProcessIsMaster();
            embeddedCrack_master.InitializeGeometry(crackMouth, crackTip);
        }

        public void InitializeGeometry(PolyLine2D initialCrack)
        {
            procs.CheckProcessIsMaster();
            embeddedCrack_master.InitializeGeometry(initialCrack);
        }

        public void Propagate(Dictionary<int, Vector> totalFreeDisplacements)
        {
            procs.CheckProcessIsMaster();
            embeddedCrack_master.Propagate(totalFreeDisplacements);
        }

        public void ScatterCrackData(IXModelMpi model)
        {
            ScatterLevelSetData(model);
            ScatterEnrichmentData(model);
        }

        public double SignedDistanceOf(XNode node) => LevelSets.SignedDistanceOf(node);

        public double SignedDistanceOf(NaturalPoint point, XContinuumElement2D element, EvalInterpolation2D interpolation)
            => LevelSets.SignedDistanceOf(point, element, interpolation);

        public void UpdateEnrichments()
        {
            procs.CheckProcessIsMaster();
            embeddedCrack_master.UpdateEnrichments();
        }

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            procs.CheckProcessIsMaster();
            embeddedCrack_master.UpdateGeometry(localGrowthAngle, growthLength);
        }

        private static (SingleCrackLsm, TipCoordinateSystem) DeserializeLevelSetData(double[] levelSetData, int subdomainID, 
            IXModelMpi model)
        {
            // Crack tip
            CartesianPoint crackTip = new CartesianPoint(levelSetData[0], levelSetData[1], levelSetData[2]);
            var tipSystem = new TipCoordinateSystem(crackTip, levelSetData[3]);

            // Nodal level set data
            XSubdomain subdomain = model.GetXSubdomain(subdomainID);
            var levelSetsBody = new Dictionary<XNode, double>();
            var levelSetsTip = new Dictionary<XNode, double>();
            int flatIdx = 4;
            foreach (XNode node in subdomain.Nodes.Values)
            {
                levelSetsBody[node] = levelSetData[flatIdx];
                levelSetsTip[node] = levelSetData[flatIdx + 1];
                flatIdx += 2;
            }
            var lsm = new SingleCrackLsm(levelSetsBody, levelSetsTip, crackTip);

            return (lsm, tipSystem);
        }

        private static void EnrichNodesAndElements(int[] bodyNodes, int[] tipNodes, int[] bodyElements, int[] tipElements,
            int subdomainID, IXModelMpi model, CrackBodyEnrichment2D bodyEnrichment, CrackTipEnrichments2D tipEnrichments)
        {
            // Clear existing enrichments
            XSubdomain subdomain = model.GetXSubdomain(subdomainID);
            foreach (XNode node in subdomain.Nodes.Values) node.EnrichmentItems.Clear();
            foreach (IXFiniteElement element in subdomain.Elements.Values) element.EnrichmentItems.Clear();

            // Enrich the nodes and elements marked by master
            foreach (int n in bodyNodes)
            {
                XNode node = subdomain.Nodes[n];
                node.EnrichmentItems[bodyEnrichment] = bodyEnrichment.EvaluateFunctionsAt(node);
            }
            foreach (int n in tipNodes)
            {
                XNode node = subdomain.Nodes[n];
                node.EnrichmentItems[tipEnrichments] = tipEnrichments.EvaluateFunctionsAt(node);
            }
            foreach (int e in bodyElements) subdomain.Elements[e].EnrichmentItems.Add(bodyEnrichment);
            foreach (int e in tipElements) subdomain.Elements[e].EnrichmentItems.Add(tipEnrichments);
        }

        private static (int[] bodyNodes, int[] tipNodes, int[] bodyElements, int[] tipElements) SerializeEnrichmentData_master(
            int subdomainID, IXModelMpi model, CrackBodyEnrichment2D bodyEnrichment, CrackTipEnrichments2D tipEnrichments)
        {
            XSubdomain subdomain = model.GetXSubdomain(subdomainID);

            // Find enriched node indices
            var bodyNodes = new List<int>();
            var tipNodes = new List<int>();
            foreach (XNode node in subdomain.Nodes.Values)
            {
                if (node.EnrichmentItems.ContainsKey(bodyEnrichment)) bodyNodes.Add(node.ID);
                if (node.EnrichmentItems.ContainsKey(tipEnrichments)) tipNodes.Add(node.ID);
            }

            // Find enriched element indices
            //TODO: This should not be necessary. Elements should not store the enrichments that affect them and instead 
            //      use the enrichment data of nodes.
            var bodyElements = new List<int>();
            var tipElements = new List<int>();
            foreach (IXFiniteElement element in subdomain.Elements.Values)
            {
                if (element.EnrichmentItems.Contains(bodyEnrichment)) bodyElements.Add(element.ID);
                if (element.EnrichmentItems.Contains(tipEnrichments)) tipElements.Add(element.ID);
            }

            return (bodyNodes.ToArray(), tipNodes.ToArray(), bodyElements.ToArray(), tipElements.ToArray());
        }

        private static double[] SerializeLevelSetData_master(SingleCrackLsm levelSets, TipCoordinateSystem tipSystem, 
            int subdomainID, IXModelMpi model)
        {
            // Serialize level set data
            XSubdomain subdomain = model.GetXSubdomain(subdomainID);
            var levelSetData = new double[2 * subdomain.NumNodes + 4];

            // Crack tip
            levelSetData[0] = levelSets.CrackTip.X;
            levelSetData[1] = levelSets.CrackTip.Y;
            levelSetData[2] = levelSets.CrackTip.Z;
            levelSetData[3] = tipSystem.RotationAngle;

            // Nodal level sets. Assume the order of enumerating them is always the same.
            int flatIdx = 4;
            foreach (XNode node in subdomain.Nodes.Values)
            {
                levelSetData[flatIdx] = levelSets.LevelSetsBody[node];
                levelSetData[flatIdx + 1] = levelSets.LevelSetsTip[node];
                flatIdx += 2;
            }

            return levelSetData;
        }

        private void ScatterEnrichmentData(IXModelMpi model)
        {
            if (procs.IsMasterProcess)
            {
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) continue;
                    int s = procs.GetSubdomainIDsOfProcess(p).First();
                    (int[] bodyNodes, int[] tipNodes, int[] bodyElements, int[] tipElements) = SerializeEnrichmentData_master(
                        s, model, CrackBodyEnrichment, CrackTipEnrichments);

                    // Send the enriched node and element indices to the corresponding process
                    MpiUtilities.SendArray<int>(procs.Communicator, bodyNodes, p, enrichmentDataTag);
                    MpiUtilities.SendArray<int>(procs.Communicator, tipNodes, p, enrichmentDataTag);
                    MpiUtilities.SendArray<int>(procs.Communicator, bodyElements, p, enrichmentDataTag);
                    MpiUtilities.SendArray<int>(procs.Communicator, tipElements, p, enrichmentDataTag);
                }
            }
            else
            {
                // Receive which nodes and elements are enriched
                int[] bodyNodes = MpiUtilities.ReceiveArray<int>(procs.Communicator, procs.MasterProcess, enrichmentDataTag);
                int[] tipNodes = MpiUtilities.ReceiveArray<int>(procs.Communicator, procs.MasterProcess, enrichmentDataTag);
                int[] bodyElements = MpiUtilities.ReceiveArray<int>(procs.Communicator, procs.MasterProcess, enrichmentDataTag);
                int[] tipElements = MpiUtilities.ReceiveArray<int>(procs.Communicator, procs.MasterProcess, enrichmentDataTag);

                int[] subdomainIds = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                EnrichNodesAndElements(bodyNodes, tipNodes, bodyElements, tipElements, subdomainIds.First(), model,
                    CrackBodyEnrichment, CrackTipEnrichments);
            }
        }

        private void ScatterLevelSetData(IXModelMpi model)
        {
            if (procs.IsMasterProcess)
            {
                this.LevelSets = embeddedCrack_master.LevelSets;
                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) continue;
                    int s = procs.GetSubdomainIDsOfProcess(p).First();
                    double[] levelSetData = SerializeLevelSetData_master(LevelSets, CrackTipEnrichments.TipSystem,
                        s, model);
                    MpiUtilities.SendArray<double>(procs.Communicator, levelSetData, p, levelSetDataTag);
                }
            }
            else
            {
                int[] subdomainIds = procs.GetSubdomainIDsOfProcess(procs.OwnRank);
                double[] levelSetData =
                    MpiUtilities.ReceiveArray<double>(procs.Communicator, procs.MasterProcess, levelSetDataTag);
                (SingleCrackLsm lsm, TipCoordinateSystem tipSystem) =
                    DeserializeLevelSetData(levelSetData, subdomainIds.First(), model);
                LevelSets = lsm;
                CrackTipEnrichments.TipSystem = tipSystem;
            }
        }
    }
}
