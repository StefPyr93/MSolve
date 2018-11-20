﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Ideally all the relevant methods should return Vector (or at least Vector for element level) 
//TODO: The map element to subdomain should be a BiList, so that it can be read in both directions and passed to native code 
//      (e.g. CUDA). It should also be cached for each element, unless there isn't enough memory or if the dof ordering is 
//      updated frequently. It can also be used for mapping vectors (e.g. ExtractVectorElementFromSubdomain), not only the
//      stiffness matrix.
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public interface ISubdomainFreeDofOrdering
    {
        DofTable FreeDofs { get; }

        int NumFreeDofs { get; }

        //TODO: This method belongs to a global vector assembler. There, it should use the element-subdomain dof map from IDofOrdering
        //TODO: What should it contain for constrained dofs?
        void AddVectorElementToSubdomain(IElement element, IVectorView elementVector, IVector subdomainVector);

        //TODO: This method belongs to a global vector assembler. There, it should use the element-subdomain dof map from IDofOrdering
        //TODO: What should it contain for constrained dofs?
        //TODOMaria: here is where the element displacements are assigned to zero if they are restrained
        double[] ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector); 

        (int[] elementDofIndices, int[] subdomainDofIndices) MapFreeDofsElementToSubdomain(IElement element);
    }
}
