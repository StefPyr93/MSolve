using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests.Tranfer
{
    internal static class TransferrerTestsData
    {
        internal static bool CheckEquality<T>(T[] array1, T[] array2)
        {
            if (array1.Length != array2.Length) return false;
            for (int i = 0; i < array1.Length; ++i)
            {
                if (!array1[i].Equals(array2[i])) return false;
            }
            return true;
        }

        internal static double[] GetArrayDataOfSubdomain(int subdomainID)
        {
            int length = subdomainID + 4;
            var data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = (subdomainID + 1) * i / 10;
            return data;
        }

        internal static SampleClass GetClassDataOfSubdomain(int subdomainID)
        {
            int length = 2 * subdomainID + 1;
            var data = new int[length];
            for (int i = 0; i < length; ++i) data[i] = subdomainID + 11;
            return new SampleClass(subdomainID, data);
        }

        internal static long GetPrimitiveDataOfSubdomain(int subdomainID) => subdomainID * 5 + 3;

        [Serializable]
        internal class SampleClass
        {
            private readonly int id;
            private readonly int[] data;
            
            public SampleClass(int id, int[] data)
            {
                this.id = id;
                this.data = data;
            }

            public int ID => id;
            public int[] Data => data;

            public int PackedArrayLength => data.Length + 1;

            public static SampleClass UnpackFromArray(int subdomainID, int[] buffer, int start, int end)
            {
                //int length = 2 * subdomainID + 1; ; // This must match GetClassDataOfSubdomain()
                int packedLength = end - start;
                var rawData = new int[packedLength - 1];
                Array.Copy(buffer, start + 1, rawData, 0, rawData.Length);
                return new SampleClass(buffer[start], rawData);
            }

            public bool Equals(SampleClass other)
            {
                if (this.id != other.id) return false;
                return CheckEquality(this.data, other.data);
            }

            public void PackIntoArray(int[] buffer, int offset)
            {
                buffer[offset] = ID;
                Array.Copy(data, 0, buffer, offset + 1, data.Length);
            }

            public override string ToString()
            {
                var builder = new StringBuilder($"id = {id}, values =");
                foreach (int d in data) builder.Append(" " + d);
                return builder.ToString();
            }
        }

        [Serializable]
        internal class SampleDto
        {
            private int[] data;

            public SampleDto(SampleClass obj)
            {
                this.data = new int[obj.Data.Length + 1];
                this.data[0] = obj.ID;
                Array.Copy(obj.Data, 0, this.data, 1, obj.Data.Length);
            }

            public SampleClass Unpack()
            {
                var rawData = new int[this.data.Length - 1];
                Array.Copy(this.data, 1, rawData, 0, rawData.Length);
                return new SampleClass(this.data[0], rawData);
            }
        }
    }
}
