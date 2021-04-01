using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

//TODO: It is possible that some previous corner nodes become internal due to the TipAdaptivePartitioner. How to handle this?
namespace ISAAR.MSolve.XFEM.Solvers
{
    public abstract class CrackedFetiDPCornerNodesBase : ICornerNodeSelection
    {
        protected readonly ICrackDescription crack;
        //private readonly Dictionary<int, HashSet<INode>> currentCornerNodes;
        protected readonly Func<ISubdomain, HashSet<INode>> getInitialCornerNodes;

        protected Dictionary<ISubdomain, bool> areCornerNodesModified;
        protected bool areGlobalCornerNodesModified = true;
        protected HashSet<INode> cornerNodesGlobal;
        protected Dictionary<ISubdomain, HashSet<INode>> cornerNodesOfSubdomains = new Dictionary<ISubdomain, HashSet<INode>>();
        protected bool isFirstAnalysis; //TODO: This should be passed to Update(). Update should not be called by the solver.

        protected CrackedFetiDPCornerNodesBase(ICrackDescription crack, Func<ISubdomain, HashSet<INode>> getInitialCornerNodes)
        {
            this.crack = crack;
            this.getInitialCornerNodes = getInitialCornerNodes;
            this.isFirstAnalysis = true;
        }

        public abstract HashSet<INode> GlobalCornerNodes { get; }

        public abstract HashSet<INode> GetCornerNodesOfSubdomain(ISubdomain subdomain);

        public abstract void Update();

        public void WriteCornerNodes(ISubdomain subdomain, string path, bool append)
        {
            using (var writer = new System.IO.StreamWriter(path, append))
            {
                int[] cornerNodes = cornerNodesOfSubdomains[subdomain].Select(n => n.ID).OrderBy(n => n).ToArray();
                writer.Write($"Corner nodes of subdomain {subdomain.ID}: ");
                foreach (int n in cornerNodes) writer.Write(n + " ");
                writer.WriteLine();
            }
        }

        protected void GatherGlobalCornerNodes()
        {
            // Define them
            cornerNodesGlobal = new HashSet<INode>();
            foreach (IEnumerable<INode> subdomainNodes in cornerNodesOfSubdomains.Values)
            {
                cornerNodesGlobal.UnionWith(subdomainNodes);
            }

            // Determine if they are modified
            areGlobalCornerNodesModified = false;
            foreach (bool modifiedSubdomain in areCornerNodesModified.Values)
            {
                if (modifiedSubdomain)
                {
                    areGlobalCornerNodesModified = true;
                    break;
                }
            }

            CornerNodeUtilities.CheckGlobalCornersAreSubdomainCorners(cornerNodesGlobal, s => cornerNodesOfSubdomains[s]);
        }

        protected void UpdateSubdomainsCorners(IEnumerable<ISubdomain> subdomains)
        {
            if (isFirstAnalysis)
            {
                foreach (ISubdomain subdomain in subdomains)
                {
                    cornerNodesOfSubdomains[subdomain] = getInitialCornerNodes(subdomain);
                }
            }

            // Keep track of subdomains with modified corner nodes.
            areCornerNodesModified = new Dictionary<ISubdomain, bool>();
            foreach (ISubdomain subdomain in cornerNodesOfSubdomains.Keys)
            {
                areCornerNodesModified[subdomain] = isFirstAnalysis;
            }

            // Remove the previous corner nodes that are no longer boundary.
            foreach (ISubdomain subdomain in cornerNodesOfSubdomains.Keys)
            {
                int numEntriesRemoved = cornerNodesOfSubdomains[subdomain].RemoveWhere(node => node.Multiplicity < 2);
                if (numEntriesRemoved > 0) areCornerNodesModified[subdomain] = true;
            }

            // Remove the previous corner nodes that do not belong to each subdomain
            foreach (ISubdomain subdomain in cornerNodesOfSubdomains.Keys)
            {
                var removedNodes = new HashSet<INode>();
                foreach (INode node in cornerNodesOfSubdomains[subdomain])
                {
                    try
                    {
                        subdomain.GetNode(node.ID);
                    }
                    catch (KeyNotFoundException)
                    {
                        removedNodes.Add(node);
                    }
                }
                foreach (INode node in removedNodes) cornerNodesOfSubdomains[subdomain].Remove(node);
                if (removedNodes.Count > 0) areCornerNodesModified[subdomain] = true;
            }

            // Add boundary Heaviside nodes and nodes of the tip element(s).
            HashSet<XNode> enrichedBoundaryNodes = FindNewEnrichedBoundaryNodes();
            foreach (XNode node in enrichedBoundaryNodes)
            {
                foreach (ISubdomain subdomain in node.SubdomainsDictionary.Values)
                {
                    cornerNodesOfSubdomains[subdomain].Add(node);
                    areCornerNodesModified[subdomain] = true;
                }
            }
        }

        private HashSet<XNode> FindNewEnrichedBoundaryNodes()
        {
            var enrichedBoundaryNodes = new HashSet<XNode>();
            foreach (CartesianPoint crackTip in crack.CrackTips)
            {
                foreach (XContinuumElement2D tipElement in crack.CrackTipElements[crackTip])
                {
                    foreach (XNode node in tipElement.Nodes)
                    {
                        if (node.Multiplicity > 1) enrichedBoundaryNodes.Add(node);
                    }
                }
            }
            foreach (var crackBodyNewNodes in crack.CrackBodyNodesNew)
            {
                foreach (XNode node in crackBodyNewNodes.Value)
                {
                    if (node.Multiplicity > 1) enrichedBoundaryNodes.Add(node);
                }
            }
            return enrichedBoundaryNodes;
        }
    }
}
