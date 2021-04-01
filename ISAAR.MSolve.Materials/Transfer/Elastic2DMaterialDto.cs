using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Transfer
{
    [Serializable]
    public class Elastic2DMaterialDto : IMaterialDto
    {
        public int id;
        public StressState2D stressState;
        public double poissonRatio;
        public double youngModulus;

        public Elastic2DMaterialDto(int id, StressState2D stressState, double youngModulus, double poissonRatio)
        {
            this.id = id;
            this.stressState = stressState;
            this.youngModulus = youngModulus;
            this.poissonRatio = poissonRatio;
        }

        public Elastic2DMaterialDto(ElasticMaterial2D material)
        {
            this.id = material.ID;
            this.stressState = material.StressState;
            this.youngModulus = material.YoungModulus;
            this.poissonRatio = material.PoissonRatio;
        }

        public int ID => id;

        public IFiniteElementMaterial Deserialize()
        {
            return new ElasticMaterial2D(stressState)
            {
                ID = this.id, YoungModulus = this.youngModulus, PoissonRatio = this.poissonRatio
            };
        }
    }
}
