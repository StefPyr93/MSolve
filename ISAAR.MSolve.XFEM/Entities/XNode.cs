using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees;

//TODO: It is not a good idea to inherit from Node. It creates many problems usually due to lack of covariance. It is better
//      to implement INode directly, but INode itself should be refactored to not expose non covariant data structures 
//      (IElement must too).
namespace ISAAR.MSolve.XFEM.Entities
{
    public class XNode : Node
    {
        public XNode(int id, double x, double y, int multiplicity = 1) : base(id, x, y, multiplicity)
        {
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }

        public XNode(int id, double x, double y, double z, int multiplicity = 1) : base(id, x, y, z, multiplicity)
        {
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }

        public Dictionary<IEnrichmentItem2D, double[]> EnrichmentItems { get; }
        public bool IsEnriched { get { return EnrichmentItems.Count > 0; } }

        public int EnrichedDofsCount
        {
            get
            {
                int count = 0;
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys) count += enrichment.Dofs.Count;
                return count;
            }
        }

        public IReadOnlyList<EnrichedDof> EnrichedDofs
        {
            get
            {
                var dofs = new List<EnrichedDof>();
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys) dofs.AddRange(enrichment.Dofs);
                return dofs;
            }
        }

        //TODO: Redesign these and their counterparts in Node. Connectivity should be done using the Dircetization interfaces.
        #region connectivity
        public new Dictionary<int, IXFiniteElement> ElementsDictionary { get; } = new Dictionary<int, IXFiniteElement>();

        public void BuildXSubdomainDictionary()
        {
            foreach (IXFiniteElement element in ElementsDictionary.Values)
            {
                SubdomainsDictionary[element.Subdomain.ID] = element.Subdomain;
            }
            Multiplicity = SubdomainsDictionary.Count;
        }
        #endregion
    }
}