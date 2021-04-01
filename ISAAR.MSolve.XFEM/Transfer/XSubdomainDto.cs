using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

//TODO: By having each process store the element factory and pass it to this object, we avoid the headache of transfering 
//      materials and especially integration rules (which use decorators to make matters worse). However, currntly only one
//      factory per subdomain is available, which means only 1 material and 1 type of integration
namespace ISAAR.MSolve.XFEM.Transfer
{
    [Serializable]
    public class XSubdomainDto //TODO: dofs, enrichments, level sets
    {
        public int id;
        public XNodeDto[] nodes;
        public XElementDto[] elements;
        public XNodalDisplacementDto[] nodalDisplacements;
        public XNodalLoadDto[] nodalLoads;

        public static XSubdomainDto CreateEmpty() => new XSubdomainDto();

        public static XSubdomainDto Serialize(XSubdomain subdomain, IDofSerializer dofSerializer)
        {
            var dto = new XSubdomainDto();
            dto.id = subdomain.ID;

            // Nodes
            dto.nodes = new XNodeDto[subdomain.NumNodes];
            int n = 0;
            foreach (XNode node in subdomain.Nodes.Values) dto.nodes[n++] = new XNodeDto(node);

            // Elements
            dto.elements = new XElementDto[subdomain.NumElements];
            int e = 0;
            foreach (IXFiniteElement element in subdomain.Elements.Values) dto.elements[e++] = new XElementDto(element);

            // Displacements
            var displacements = new List<XNodalDisplacementDto>();
            foreach (XNode node in subdomain.Nodes.Values)
            {
                foreach (Constraint constraint in node.Constraints)
                {
                    displacements.Add(new XNodalDisplacementDto(node, constraint, dofSerializer));
                }
            }
            dto.nodalDisplacements = displacements.ToArray();

            // Nodal loads
            dto.nodalLoads = new XNodalLoadDto[subdomain.NodalLoads.Count];
            for (int i = 0; i < subdomain.NodalLoads.Count; ++i)
            {
                dto.nodalLoads[i] = new XNodalLoadDto(subdomain.NodalLoads[i], dofSerializer);
            }

            return dto;
        }

        public XSubdomain Deserialize(IDofSerializer dofSerializer, IXFiniteElementFactory elementFactory)
        {
            var subdomain = new XSubdomain(id);

            // Nodes
            foreach (XNodeDto n in this.nodes)
            {
                XNode node = n.Deserialize();
                subdomain.Nodes[node.ID] = node;
            }

            // Elements
            foreach (XElementDto e in this.elements)
            {
                IXFiniteElement element = e.Deserialize(elementFactory, id => subdomain.Nodes[id]);
                subdomain.Elements.Add(element.ID, element);
            }

            // Displacements
            foreach (XNodalDisplacementDto d in this.nodalDisplacements) d.Deserialize(subdomain.Nodes, dofSerializer);

            // Nodal loads
            foreach (XNodalLoadDto nl in this.nodalLoads) subdomain.NodalLoads.Add(nl.Deserialize(subdomain.Nodes, dofSerializer));

            return subdomain;
        }
    }
}
