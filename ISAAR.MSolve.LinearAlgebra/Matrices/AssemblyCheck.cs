using ISAAR.MSolve.LinearAlgebra.Commons;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public static class AssemblyCheck
    {
        public static bool WritingPhase { get; set; } = false;

        public static bool isSecondAssembly { get; set; } = false;

        //public static bool isSubdomainId { get; set; }

        public static Dictionary<int,double[,]> elementMatrices = new Dictionary<int,double[,]>();

        public static Dictionary<int, bool> areMatricesSame = new Dictionary<int, bool>();

        public static void ProcessElementStiffnessMatrix(int elementId, double[,] matrix)
        {
            if(WritingPhase)
            {
                if(isSecondAssembly)
                {
                    elementMatrices.Add(elementId,matrix);
                }
            }
            else
            {
                if (isSecondAssembly)
                {
                    areMatricesSame[elementId] = AreDisplacementsSame(elementMatrices[elementId], matrix);
                }
            }
        }

        public static bool AreDisplacementsSame(double[,] expectedValues,
            double[,] computedValues, double tol =1E-13)
        {
            var comparer = new ValueComparer(tol);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public static void ViewResults()
        {
            var a = areMatricesSame;
        }


    }
}
