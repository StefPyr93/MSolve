using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Entities
{
    public class Cluster
    {
        public Cluster(int id)
        {
            this.ID = id;
            this.Subdomains = new List<ISubdomain>();
        }

        public int ID { get; }

        public IList<ISubdomain> Subdomains { get; }
    }
}
