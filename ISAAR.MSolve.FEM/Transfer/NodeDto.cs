using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Transfer
{
    [Serializable]
    public class NodeDto
    {
        public int id;

        /// <summary>
        /// Transfering all the subdomain IDs is both time consuming and will cause problems since not all subdomains are 
        /// available to each process.
        /// </summary>
        public int multiplicity; //TODO: Storing and working with a list of boundary nodes of each subdomain would be cleaner.

        public double x, y, z;

        public NodeDto(int id, double x, double y, double z, int multiplicity)
        {
            this.id = id;
            this.multiplicity = multiplicity;
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public NodeDto(Node node)
        {
            this.id = node.ID;
            this.multiplicity = node.Multiplicity;
            this.x = node.X;
            this.y = node.Y;
            this.z = node.Z;
        }

        //public static NodeDto Serialize(Node node)
        //{
        //    var trans = new NodeDto();
        //    trans.id = node.ID;
        //    trans.x = node.X;
        //    trans.y = node.Y;
        //    trans.z = node.Z;
        //    return trans;
        //}

        public Node Deserialize() => new Node(id, x, y, z, multiplicity);
    }
}
