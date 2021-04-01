using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Entities;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Distributed;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ModelMpiRedundant : ModelMpiRedundantBase<Model>
    {
        public ModelMpiRedundant(ProcessDistribution processDistribution, Func<Model> createModel) 
            : base(processDistribution, createModel)
        {
        }
    }
}
