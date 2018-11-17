﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IStructuralModel_v2
    {
		Dictionary<int, ISubdomain_v2> ISubdomainsDictionary { get;  }
	    IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; }
		void AssignLoads();
		void AssignMassAccelerationHistoryLoads(int timeStep);
	}
}
