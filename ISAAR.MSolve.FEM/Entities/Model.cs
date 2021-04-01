﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Transfer;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: find what is going on with the dynamic loads and refactor them. That 564000000 in AssignMassAccelerationHistoryLoads()
//      cannot be correct.
//TODO: ConnectDataStructures() should not be called twice. There should be a flag that determines if it has been called. If it
//      has, the method should just return without doing anything.
//TODO: Replace all IList with IReadOnlyList. Even better, have a different class to create the model than the one used to 
//      store the entities, so that they can be accessed by solvers, analyzers & loggers. Could the latter be the same for FEM, 
//      IGA, XFEM?
namespace ISAAR.MSolve.FEM.Entities
{
    public class Model : IModel
    {
        //public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

        public Dictionary<int, Cluster> Clusters { get; } = new Dictionary<int, Cluster>();

        public Table<INode, IDofType, double> Constraints { get; private set; } = new Table<INode, IDofType, double>();//TODOMaria: maybe it's useless in model class

        public IDofSerializer DofSerializer { get; set; } = new StandardDofSerializer();

        public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

        public IList<ElementMassAccelerationHistoryLoad> ElementMassAccelerationHistoryLoads { get; } 
            = new List<ElementMassAccelerationHistoryLoad>();
        public IList<ElementMassAccelerationLoad> ElementMassAccelerationLoads { get; } 
            = new List<ElementMassAccelerationLoad>();
        public IList<Load> Loads { get; private set; } = new List<Load>();
        public IList<MassAccelerationLoad> MassAccelerationLoads { get; } = new List<MassAccelerationLoad>();
        public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } = new List<IMassAccelerationHistoryLoad>();

        public Dictionary<int, Node> NodesDictionary { get; } = new Dictionary<int, Node>();

        public int NumElements => ElementsDictionary.Count;
        public int NumNodes => NodesDictionary.Count;
        public int NumSubdomains => SubdomainsDictionary.Count;

        public Dictionary<int, Subdomain> SubdomainsDictionary { get; } = new Dictionary<int, Subdomain>();

        public IList<ITimeDependentNodalLoad> TimeDependentNodalLoads { get; private set; } = new List<ITimeDependentNodalLoad>();

        public IGlobalFreeDofOrdering GlobalDofOrdering { get; set; }

        //TODO: nodal loads and prescribed displacements are calculated by dedicated objects. The rest should follow asap.
        public void ApplyLoads()
        {
            foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.Forces.Clear();
            AssignElementMassLoads();
            AssignMassAccelerationLoads();
        }

        public void ApplyMassAccelerationHistoryLoads(int timeStep)
        {
            if (MassAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(MassAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in MassAccelerationHistoryLoads)
                {
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });
                }

                foreach (Subdomain subdomain in SubdomainsDictionary.Values)
                {
                    foreach (Element element in subdomain.Elements.Values)
                    {
                        double[] accelerationForces = element.ElementType.CalculateAccelerationForces(element, m);
                        subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, accelerationForces, subdomain.Forces);
                    }
                }
            }

            foreach (ElementMassAccelerationHistoryLoad load in ElementMassAccelerationHistoryLoads)
            {
                MassAccelerationLoad hl = new MassAccelerationLoad()
                {
                    Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF
                };
                Element element = load.Element;
                ISubdomain subdomain = element.Subdomain;
                var accelerationForces = element.ElementType.CalculateAccelerationForces(
                    load.Element, (new MassAccelerationLoad[] { hl }).ToList());
                GlobalDofOrdering.GetSubdomainDofOrdering(subdomain).AddVectorElementToSubdomain(element, accelerationForces,
                    subdomain.Forces);
            }
        }

        //TODO: These are not distributed correctly. They should be handled like steady state nodal loads.
        public void AssignTimeDependentNodalLoads(int timeStep)
        {
            foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
            {
                double amountPerSubdomain = load.GetLoadAmount(timeStep) / load.Node.Multiplicity;
                foreach (ISubdomain subdomain in load.Node.SubdomainsDictionary.Values)
                {
                    int idx = subdomain.FreeDofOrdering.FreeDofs[load.Node, load.DOF];
                    subdomain.Forces[idx] += amountPerSubdomain;
                }
            }
        }

        //What is the purpose of this method? If someone wanted to clear the Model, they could just create a new one.
        public void Clear()
        {
            Loads.Clear();
            Clusters.Clear();
            SubdomainsDictionary.Clear();
            ElementsDictionary.Clear();
            NodesDictionary.Clear();
            GlobalDofOrdering = null;
            Constraints.Clear();
            ElementMassAccelerationHistoryLoads.Clear();
            ElementMassAccelerationLoads.Clear();
            MassAccelerationHistoryLoads.Clear();
            MassAccelerationLoads.Clear();
        }

        // Warning: This is called by the analyzer, so that the user does not have to call it explicitly. However, it is must be 
        // called explicitly before the AutomaticDomainDecompositioner is used.
        public void ConnectDataStructures()
        {         
            BuildInterconnectionData();
            AssignConstraints();
            AssignNodalLoadsToSubdomains();
            RemoveInactiveTransientNodalLoads();

            //TODOSerafeim: This should be called by the analyzer, which defines when the dofs are ordered and when the global vectors/matrices are built.
            //AssignLoads();
        }

        public IEnumerable<IElement> EnumerateElements() => ElementsDictionary.Values;
        public IEnumerable<INode> EnumerateNodes() => NodesDictionary.Values;
        public IEnumerable<ISubdomain> EnumerateSubdomains() => SubdomainsDictionary.Values;

        public IElement GetElement(int elementID) => ElementsDictionary[elementID];
        public INode GetNode(int nodeID) => NodesDictionary[nodeID];
        public ISubdomain GetSubdomain(int subdomainID) => SubdomainsDictionary[subdomainID];

        //TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
        //      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
        //      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
        //      It is too easy to access the wrong instance of the constraint. 
        private void AssignConstraints()
        {
            foreach (Node node in NodesDictionary.Values)
            {
                if (node.Constraints == null) continue;
                foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            }

            foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.ExtractConstraintsFromGlobal(Constraints);
        }

        private void AssignElementMassLoads()
        {
            foreach (ElementMassAccelerationLoad load in ElementMassAccelerationLoads)
            {
                ISubdomain subdomain = load.Element.Subdomain;
                var accelerationForces = load.Element.ElementType.CalculateAccelerationForces(
                    load.Element, MassAccelerationLoads);
                GlobalDofOrdering.GetSubdomainDofOrdering(subdomain).AddVectorElementToSubdomain(load.Element,
                    accelerationForces, subdomain.Forces);
            }
        }

        private void AssignMassAccelerationLoads()
        {
            if (MassAccelerationLoads.Count < 1) return;

            foreach (Subdomain subdomain in SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.Elements.Values)
                {
                    subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element,
                        element.ElementType.CalculateAccelerationForces(element, MassAccelerationLoads),
                        subdomain.Forces);
                }
            }
        }

        private void AssignNodalLoadsToSubdomains()
        {
            // Remove inactive loads added by the user
            var activeLoads = new List<Load>(Loads.Count);
            foreach (Load load in Loads)
            {
                bool isConstrained = Constraints.Contains(load.Node, load.DOF);
                if (!isConstrained) activeLoads.Add(load);
            }
            Loads = activeLoads;

            // Assign the rest to their subdomains without scaling them. That will be done later by the analyzer and solver.
            foreach (Load load in Loads)
            {
                foreach (Subdomain subdomain in load.Node.SubdomainsDictionary.Values) subdomain.NodalLoads.Add(load);
            }
        }

        private void BuildElementDictionaryOfEachNode()
        {
            foreach (Element element in ElementsDictionary.Values)
            {
                foreach (Node node in element.Nodes) node.ElementsDictionary[element.ID] = element;
            }
        }

        private void BuildInterconnectionData()//TODOMaria: maybe I have to generate the constraints dictionary for each subdomain here
        {
            BuildSubdomainOfEachElement();
            //DuplicateInterSubdomainEmbeddedElements();
            BuildElementDictionaryOfEachNode();
            foreach (Node node in NodesDictionary.Values) node.BuildSubdomainDictionary();

            //BuildNonConformingNodes();

            foreach (Subdomain subdomain in SubdomainsDictionary.Values) subdomain.DefineNodesFromElements();
        }

        private void BuildSubdomainOfEachElement()
        {
            foreach (Subdomain subdomain in SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.Elements.Values) element.Subdomain = subdomain;
            }
        }

        private void BuildNonConformingNodes()
        {
            List<int> subIDs = new List<int>();
            foreach (Element element in ElementsDictionary.Values)
            {
                subIDs.Clear();

                foreach (Node node in element.Nodes)
                {
                    foreach (int subID in node.SubdomainsDictionary.Keys)
                    {
                        if (!subIDs.Contains(subID)) subIDs.Add(subID);

                    }
                }

                foreach (Node node in element.Nodes)
                {
                    foreach (int subID in subIDs)
                    {
                        if (!node.SubdomainsDictionary.ContainsKey(subID))
                        {
                            node.NonMatchingSubdomainsDictionary.Add(subID, SubdomainsDictionary[subID]);
                        }
                    }
                }
                
            }
        }

        private void DuplicateInterSubdomainEmbeddedElements()
        {
            foreach (var e in ElementsDictionary.Values.Where(x => x.ElementType is IEmbeddedElement))
            {
                var subs = ((IEmbeddedElement)e.ElementType).EmbeddedNodes.Select(x => x.EmbeddedInElement.Subdomain).Distinct();
                foreach (var s in subs.Where(x => x.ID != e.Subdomain.ID)) s.Elements.Add(e.ID, e);
            }
        }

        private void RemoveInactiveTransientNodalLoads()
        {
            var activeLoadsDynamic = new List<ITimeDependentNodalLoad>(TimeDependentNodalLoads.Count);
            foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
            {
                bool isConstrained = Constraints.Contains(load.Node, load.DOF);
                if (!isConstrained) activeLoadsDynamic.Add(load);
            }
            TimeDependentNodalLoads = activeLoadsDynamic;
        }

        //private void EnumerateGlobalDOFs()
        //{
        //    totalDOFs = 0;
        //    Dictionary<int, List<DOFType>> nodalDOFTypesDictionary = new Dictionary<int, List<DOFType>>();
        //    foreach (Element element in elementsDictionary.Values)
        //    {
        //        for (int i = 0; i < element.Nodes.Count; i++)
        //        {
        //            if (!nodalDOFTypesDictionary.ContainsKey(element.Nodes[i].ID))
        //                nodalDOFTypesDictionary.Add(element.Nodes[i].ID, new List<DOFType>());
        //            nodalDOFTypesDictionary[element.Nodes[i].ID].AddRange(element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
        //        }
        //    }

        //    foreach (Node node in nodesDictionary.Values)
        //    {
        //        Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
        //        foreach (DOFType dofType in nodalDOFTypesDictionary[node.ID].Distinct())
        //        {
        //            int dofID = 0;
        //            #region removeMaria
        //            //foreach (DOFType constraint in node.Constraints)
        //            //{
        //            //    if (constraint == dofType)
        //            //    {
        //            //        dofID = -1;
        //            //        break;
        //            //    }
        //            //}
        //            #endregion

        //            foreach (var constraint in node.Constraints)
        //            {
        //                if (constraint.DOF == dofType)
        //                {
        //                    dofID = -1;
        //                    break;
        //                }
        //            }

        //            //// TODO: this is not applicable! Embedded nodes have to do with the embedded element and not with the host
        //            //// User should define which DOFs would be dependent on the host element. For our case
        //            //// we should select between translational and translational + rotational
        //            //var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
        //            ////if (node.EmbeddedInElement != null && node.EmbeddedInElement.ElementType.GetDOFTypes(null)
        //            ////    .SelectMany(d => d).Count(d => d == dofType) > 0)
        //            ////    dofID = -1;
        //            //if (embeddedNode != null && embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null)
        //            //    .SelectMany(d => d).Count(d => d == dofType) > 0)
        //            //    dofID = -1;

        //            if (dofID == 0)
        //            {
        //                dofID = totalDOFs;
        //                totalDOFs++;
        //            }
        //            dofsDictionary.Add(dofType, dofID);
        //        }

        //        nodalDOFsDictionary.Add(node.ID, dofsDictionary);
        //    }
        //}
    }
}
