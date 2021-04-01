using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Transfer
{
    public interface IMaterialDto
    {
        int ID { get; }
        IFiniteElementMaterial Deserialize();
    }
}
