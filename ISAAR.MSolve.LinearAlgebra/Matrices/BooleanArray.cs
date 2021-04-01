using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public class BooleanArray
    {
        bool[] booleanArray;

        public BooleanArray(int size)
        {
            booleanArray = new bool[size];
        }

        public void SetTrue(int position)
        {
            booleanArray[position] = true;
        }

        public bool isTrue(int position)
        {
            return booleanArray[position];
        }
    }
}
