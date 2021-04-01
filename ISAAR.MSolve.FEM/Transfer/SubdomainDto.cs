using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Transfer.Elements;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Materials.Transfer;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Transfer
{
    [Serializable]
    public class SubdomainDto
    {
        public int id;
        public NodeDto[] nodes;
        public IElementDto[] elements;
        public IMaterialDto[] materials;
        public NodalDisplacementDto[] nodalDisplacements;
        public NodalLoadDto[] nodalLoads;

        public static SubdomainDto CreateEmpty() => new SubdomainDto();

        public static SubdomainDto Serialize(Subdomain subdomain, IDofSerializer dofSerializer)
        {
            var dto = new SubdomainDto();
            dto.id = subdomain.ID;

            // Nodes
            dto.nodes = new NodeDto[subdomain.NumNodes];
            int n = 0;
            foreach (Node node in subdomain.Nodes.Values) dto.nodes[n++] = new NodeDto(node);

            // Elements
            dto.elements = new IElementDto[subdomain.NumElements];
            var elementSerializer = new ElementSerializer();
            int e = 0;
            foreach (Element element in subdomain.Elements.Values) dto.elements[e++] = elementSerializer.Serialize(element);

            // Materials
            // More than 1 elements may have the same material properties. First gather the unique ones.
            var uniqueMaterials = new Dictionary<int, IFiniteElementMaterial>();
            foreach (Element element in subdomain.Elements.Values)
            {
                // Each element is assumed to have the same material at all GPs.
                IFiniteElementMaterial elementMaterial = element.ElementType.Materials[0];
                uniqueMaterials[elementMaterial.ID] = elementMaterial;
            }
            dto.materials = new IMaterialDto[uniqueMaterials.Count];
            var materialSerializer = new MaterialSerializer();
            int counter = 0;
            foreach (IFiniteElementMaterial material in uniqueMaterials.Values)
            {
                dto.materials[counter++] = materialSerializer.Serialize(material);
            }

            // Displacements
            var displacements = new List<NodalDisplacementDto>();
            foreach (Node node in subdomain.Nodes.Values)
            {
                foreach (Constraint constraint in node.Constraints)
                {
                    displacements.Add(new NodalDisplacementDto(node, constraint, dofSerializer));
                }
            }
            dto.nodalDisplacements = displacements.ToArray();

            // Nodal loads
            dto.nodalLoads = new NodalLoadDto[subdomain.NodalLoads.Count];
            for (int i = 0; i < subdomain.NodalLoads.Count; ++i)
            {
                dto.nodalLoads[i] = new NodalLoadDto(subdomain.NodalLoads[i], dofSerializer);
            }

            return dto;
        }

        public Subdomain Deserialize(IDofSerializer dofSerializer)
        {
            var subdomain = new Subdomain(this.id);

            // Nodes
            foreach (NodeDto n in this.nodes)
            {
                Node node = n.Deserialize();
                subdomain.Nodes[node.ID] = node;
            }

            // Materials
            var allMaterials = new Dictionary<int, IFiniteElementMaterial>();
            foreach (IMaterialDto m in this.materials) allMaterials[m.ID] = m.Deserialize();

            // Elements
            foreach (IElementDto e in this.elements)
            {
                Element element = e.Deserialize(subdomain.Nodes, allMaterials);
                subdomain.Elements.Add(element.ID, element);
            }

            // Displacements
            foreach (NodalDisplacementDto d in this.nodalDisplacements) d.Deserialize(subdomain.Nodes, dofSerializer);

            // Nodal loads
            foreach (NodalLoadDto nl in this.nodalLoads) subdomain.NodalLoads.Add(nl.Deserialize(subdomain.Nodes, dofSerializer));

            return subdomain;
        }
    }
}
