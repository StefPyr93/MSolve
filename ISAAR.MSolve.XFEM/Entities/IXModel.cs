using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;

namespace ISAAR.MSolve.XFEM.Entities
{
    public interface IXModelMpi : IModelMpi
    {
        IDomain2DBoundary Boundary { get; }

        XModel RawModel { get; }

        XSubdomain GetXSubdomain(int subdomainID);

        void ScatterSubdomains(HashSet<int> modifiedSubdomains);
        void ScatterSubdomainsState();
    }
}
