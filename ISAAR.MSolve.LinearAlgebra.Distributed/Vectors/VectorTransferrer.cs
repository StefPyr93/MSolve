using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using MPI;

//TODO: Repeat many of the methods of ISubdomainDataTransferrer here, but only for Vector. These will call some of the methods
//      in ISubdomainDataTransferrer.
//TODO: Perhaps I should provide overloads for the summation method that accept a delegate that calculates the subdomain vectors 
//      instead of the vectors themselves. 
//TODO: Also provide overloads to afacilitate the operation y = A * x + y . It is important here to avoid extra allocations and 
//      additions, since the subdomain vectors usually come from a multiplication like this. E.g. define a GEMV method that 
//      takes adds the multiplication result on top of the rhs vector, which is also passed in 
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Vectors
{
    public class VectorTransferrer
    {

        private readonly ProcessDistribution procs;

        public VectorTransferrer(ProcessDistribution processDistribution)
        {
            this.procs = processDistribution;
        }

        /// <summary>
        /// If the length of the vector is already known in each process, use <see cref="BroadcastVector(ref Vector, int)"/> 
        /// instead.
        /// </summary>
        /// <param name="vector"></param>
        public void BroadcastVector(ref Vector vector)
        {
            int length = -1;
            if (procs.IsMasterProcess) length = vector.Length;
            procs.Communicator.Broadcast<int>(ref length, procs.MasterProcess);
            BroadcastVector(ref vector, length);
        }

        /// <summary>
        /// More efficient than <see cref="BroadcastVector(ref Vector)"/>, but each process must know the length of the vector
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="length"></param>
        public void BroadcastVector(ref Vector vector, int length)
        {
            double[] data;
            if (procs.IsMasterProcess) data = vector.RawData;
            else data = new double[length];
            procs.Communicator.Broadcast<double>(ref data, procs.MasterProcess);
            if (!procs.IsMasterProcess) vector = Vector.CreateFromArray(data);
        }

        /// <summary>
        /// In master process the sum will be returned. In other processes, null will be returned.
        /// </summary>
        /// <param name="processVectors"></param>
        public Vector SumVectors(IEnumerable<Vector> processVectors)
        {
            Vector totalVector_master = null;
            if (procs.IsMasterProcess)
            {
                Vector first = processVectors.First();
                totalVector_master = Vector.CreateZero(first.Length);
            }
            SumVectors(processVectors, totalVector_master);
            return totalVector_master;
        }

        public void SumVectors(IEnumerable<Vector> processVectors, Vector totalVector_master)
        {
            // First add the vectors of this process
            Vector processSum = null;
            using (IEnumerator<Vector> iterator = processVectors.GetEnumerator())
            {
                bool isNotEmpty = iterator.MoveNext();
                if (!isNotEmpty) throw new ArgumentException("There must be at least 1 vector");
                processSum = iterator.Current.Copy(); // Optimization for the first one
                while (iterator.MoveNext()) processSum.AddIntoThis(iterator.Current); 
            }

            // MPI-reduce the vectors of all processes
            double[] asArray = processSum.RawData;
            double[] totalVectorData_master = null;
            if (procs.IsMasterProcess) totalVectorData_master = totalVector_master.RawData;
            procs.Communicator.Reduce<double>(asArray, Operation<double>.Add, procs.MasterProcess, ref totalVectorData_master);
        }
    }
}
