using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: In general I need a class to transfer continuum elements, with all posible combinations of interpolations, integrations
//      and different material per gauss point. Going further, I need a class to transfer all possible elements 
//      (continuum, structural, thermal,2 2D, 3D, etc).
namespace ISAAR.MSolve.FEM.Transfer.Elements
{
    [Serializable]
    public class Quad4Dto : IElementDto
    {
        public int id;
        public int material;
        public int node0;
        public int node1;
        public int node2;
        public int node3;
        public double thickness;

        public Quad4Dto(int id, int material, int[] nodes, double thickness)
        {
            this.id = id;
            this.material = material;
            this.node0 = nodes[0];
            this.node1 = nodes[1];
            this.node2 = nodes[2];
            this.node3 = nodes[3];
            this.thickness = thickness;
        }

        public Quad4Dto(Element element)
        {
            this.id = element.ID;
            this.node0 = element.Nodes[0].ID;
            this.node1 = element.Nodes[1].ID;
            this.node2 = element.Nodes[2].ID;
            this.node3 = element.Nodes[3].ID;

            var quad4 = ((ContinuumElement2D)element.ElementType);
            this.thickness = quad4.Thickness;
            this.material = quad4.Materials[0].ID; // Assuming that the (initial at least) properties are the same across all GPs
        }

        public Element Deserialize(IReadOnlyDictionary<int, Node> allNodes, Dictionary<int, IFiniteElementMaterial> allMaterials)
        {
            var factory = new ContinuumElement2DFactory(0, null, null);
            Node[] elemNodes = { allNodes[node0], allNodes[node1], allNodes[node2], allNodes[node3] };
            ElasticMaterial2D elemMaterial = (ElasticMaterial2D)allMaterials[material]; //TODO: Use ContinuumMaterial in continuum elements
            var elemType = factory.CreateElement(CellType.Quad4, elemNodes, thickness, elemMaterial, null); //TODO: also thisfer dynamic properties

            var elemWrapper = new Element() { ID = id, ElementType = elemType };
            foreach (Node node in elemNodes) elemWrapper.AddNode(node);
            return elemWrapper;
        }

        //public static Quad4Dto Serialize(Element element, int material, double thickness)
        //{
        //    var trans = new Quad4Dto();
        //    trans.id = element.ID;
        //    trans.material = material;
        //    trans.node0 = element.Nodes[0].ID;
        //    trans.node1 = element.Nodes[1].ID;
        //    trans.node2 = element.Nodes[2].ID;
        //    trans.node3 = element.Nodes[3].ID;
        //    trans.thickness = thickness;
        //    return trans;
        //}


        //public static Element[] Deserialize(Quad4Dto[] elements, 
        //    Dictionary<int, Node> allNodes, Dictionary<int, IContinuumMaterial3D> allMaterials)
        //{
        //    var resultElements = new Element[elements.Length];
        //    var factory = new ContinuumElement2DFactory(0, null, null);
        //    for (int e = 0; e < elements.Length; ++e)
        //    {
        //        Quad4Dto trans = elements[e];
        //        Node[] elemNodes =
        //        {
        //            allNodes[trans.node0], allNodes[trans.node1], allNodes[trans.node2], allNodes[trans.node3]
        //        };
        //        ElasticMaterial2D elemMaterial = (ElasticMaterial2D)allMaterials[trans.material]; //TODO: Use ContinuumMaterial in continuum elements
        //        var elemType = factory.CreateElement(CellType.Quad4, elemNodes, trans.thickness, elemMaterial, null); //TODO: also transfer dynamic properties

        //        var elemWrapper = new Element() { ID = trans.id, ElementType = elemType };
        //        foreach (Node node in elemNodes) elemWrapper.AddNode(node);
        //        resultElements[e] = elemWrapper;
        //    }
        //    return resultElements;
        //}
    }
}
