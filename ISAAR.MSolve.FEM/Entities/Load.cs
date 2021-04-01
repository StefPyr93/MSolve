﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Many boundary (including this) should depend on INode and IElement, in order to work for all discretization methods
namespace ISAAR.MSolve.FEM.Entities
{
    public class Load : INodalLoad
    {
        INode INodalLoad.Node => Node;
        public Node Node { get; set; }

        public IDofType DOF { get; set; }

        public double Amount { get; set; }

    }
}
