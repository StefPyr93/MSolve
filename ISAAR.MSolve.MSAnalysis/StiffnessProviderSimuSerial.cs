using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.MultiscaleAnalysis;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge
{
    /// <summary>
    /// Element stiffness matrix provider for simultaneous calculation of global stiffness matrix macroscopic variables in multiscale FE2 scheme
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class StiffnessProviderSimuSerial : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members
        private SubdomainCalculationsAndAssemblySerial host;

        public StiffnessProviderSimuSerial(SubdomainCalculationsAndAssemblySerial host)
        {
            this.host = host;
        }

        public IMatrix Matrix(IElement element) //TODOGer IMatrix2D will be changed to Matrix etc.
        {
            var elementMatrix = element.ElementType.StiffnessMatrix(element);
            host.UpdateVectors(element, elementMatrix);
            return elementMatrix;
        }

        #endregion
    }
}






    
        //#region IElementMatrixProvider Members

        //public IMatrix2D Matrix(IElement element)
        //{
        //    return element.IElementType.StiffnessMatrix(element);
        //}

        //#endregion
    
