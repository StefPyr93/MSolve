using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IFiniteElementMaterial: ICloneable
    {
        //this is only for Structural elements with strain and stress tensors.

        /// <summary>
        /// This should be used to differentiate the materials used in the same model instead of each material class 
        /// having a unique code. E.g. if 2 elastic materials with different E are used, then they get ids 0 and 1. 
        /// To differentiate material classes from each other, we can use material.GetType().
        /// </summary>
        int ID { get; } 

        bool Modified { get; }
        void ResetModified();
        double[] Coordinates { get; set; }
        double YoungModulus { get; }
        double PoissonRatio { get; }

        void SaveState();
        void ClearState();
        void ClearStresses();
    }
}
