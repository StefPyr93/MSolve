using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using MPI;

//TODO: Ideally these would be strategy objects
//TODO: Send/Receive/BroadcastArray should be in a dedicated ArrayTransferrer class, similar to VectorTransferrer and 
//      MatrixTransferrer.
namespace ISAAR.MSolve.LinearAlgebra.Distributed
{
    public static class MpiUtilities
    {
        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="comm"></param>
        /// <param name="values">
        /// In <paramref name="root"/> it must contain the array to broadcast. In other processes it can be null before calling  
        /// this method and will be allocated and filled when the method returns.
        /// </param>
        /// <param name="root"></param>
        public static void BroadcastArray<T>(Intracommunicator comm, ref T[] values, int root)
        {
            // First broadcast the length of the array
            int rank = comm.Rank;
            int length = -1;
            if (rank == root) length = values.Length;
            comm.Broadcast<int>(ref length, root);

            // Allocate buffers in receiving processes
            if (rank != root) values = new T[length];

            //The broadcast the whole array
            comm.Broadcast<T>(ref values, root);
        }

        public static void BroadcastVector(this Intracommunicator comm, ref Vector vector, int length, int root)
        {
            //TODO: Use a dedicated class for MPI communication of Vector. This class belongs to a project LinearAlgebra.MPI.
            //      Avoid copying the array.
            double[] asArray = null;
            if (comm.Rank == root) asArray = vector.CopyToArray();
            else asArray = new double[length];
            comm.Broadcast<double>(ref asArray, root);
            vector = Vector.CreateFromArray(asArray);
        }

        public static void DoInTurn(Intracommunicator comm, Action action)
        {
            comm.Barrier();
            int token = 0;
            if (comm.Rank == 0)
            {
                // Print error message
                action();

                // Send token to our right neighbor
                comm.Send(token, (comm.Rank + 1) % comm.Size, 0);

                // Receive token from our left neighbor
                comm.Receive((comm.Rank + comm.Size - 1) % comm.Size, 0, out token);
            }
            else
            {
                // Receive token from our left neighbor
                comm.Receive((comm.Rank + comm.Size - 1) % comm.Size, 0, out token);

                // Print error message
                action();

                // Pass on the token to our right neighbor
                comm.Send(token, (comm.Rank + 1) % comm.Size, 0);
            }
            comm.Barrier();
        }

        //TODO: Perhaps client code (in root) can work with vectors that have the gathered flattened array as their backing end. 
        //      I need dedicated classes for these (e.g. OffsetVector)
        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="comm"></param>
        /// <param name="arrayToGather">
        /// For root process, whatever is passed in will be returned in the result. If calculating it takes time, root can 
        /// pass in null instead.
        /// </param>
        /// <param name="arrayLengths_root"></param>
        /// <param name="root"></param>
        /// <returns></returns>
        public static T[][] GatherArrays<T>(Intracommunicator comm, T[] arrayToGather, 
            int[] arrayLengths_root, int root)
        {
            //if (comm.Rank == root)
            //{
            //    Console.Write($"Process {comm.Rank}: Counts: ");
            //    for (int p = 0; p < comm.Size; ++p) Console.Write(arrayLengths_root[p] + " ");
            //    Console.WriteLine();
            //}
            //DoInTurn(comm, () => Console.WriteLine($"Process {comm.Rank}: array.Length = {arrayToGather.Length}"));

            //TODO: perhaps do optimizations for the array of the root process
            if (comm.Rank == root) Debug.Assert(arrayLengths_root.Length == comm.Size,
                $"There are {arrayLengths_root.Length} arrays, but {comm.Size} processes."); //TODO: this will cause a deadlock in other processes

            // Gather all values to root
            T[] allValues = comm.GatherFlattened<T>(arrayToGather, arrayLengths_root, root);

            // Copy the values to individual arrays
            T[][] processArrays = null;
            if (comm.Rank == root)
            {
                processArrays = new T[comm.Size][];
                int offset = 0;
                for (int p = 0; p < comm.Size; ++p)
                {
                    if (p == root) processArrays[root] = arrayToGather;
                    int length = arrayLengths_root[p];
                    processArrays[p] = new T[length];
                    Array.Copy(allValues, offset, processArrays[p], 0, length);
                    offset += length;
                }
            }

            return processArrays;
        }

        

        //public static T[] GatherFromSubdomains<T>(Intracommunicator comm, T subdomainData, int masterProcess)
        //{
        //Dictionary<int, DofTable> subdomainCornerDofOrderings = null;
        //var tableSerializer = new DofTableSerializer(dofSerializer);
        //if (rank == masterProcess)
        //{
        //    subdomainCornerDofOrderings = new Dictionary<int, DofTable>();
        //    subdomainCornerDofOrderings[masterProcess] = SubdomainDofs.CornerDofOrdering;
        //    var requests = new Dictionary<int, ReceiveRequest>();
        //    var pending = new HashSet<int>();
        //    for (int p = 0; p < comm.Size; ++p)
        //    {
        //        if (p == masterProcess) continue;
        //        requests[p] = comm.ImmediateReceive<int[]>(p, cornerDofOrderingTag);
        //        pending.Add(p);
        //    }

        //    while (pending.Count > 0) // perhaps the thread should sleep, until a request is completed
        //    {
        //        foreach (int p in pending)
        //        {
        //            CompletedStatus status = requests[p].Test();
        //            if (status != null)
        //            {

        //            }
        //        }
        //    }

        //    for (int p = 0; p < comm.Size; ++p)
        //    {
        //        if (p == masterProcess) continue;
        //        int[] flatTable = requests[p].Get
        //    }
        //}
        //else
        //{

        //    comm.ImmediateSend(SubdomainDofs.CornerDofOrdering)
        //}
        //}

        /// <summary>
        /// Must be paired with <see cref="SendArray(Intracommunicator, int[], int, int)(Intracommunicator, int, int)"/>.
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="source"></param>
        /// <param name="tag"></param>
        /// <returns></returns>
        public static T[] ReceiveArray<T>(Intracommunicator comm, int source, int tag)
        {
            // Receiving the length is needed to allocate a buffer before receiving the whole array. Furthermore it 
            // blocks all processes except for the one currently processed master, preventing them from sending the whole 
            // arrays. Not sure if this is good or bad.
            int length = comm.Receive<int>(source, tag);
            var buffer = new T[length];
            comm.Receive<T>(source, tag, ref buffer);
            return buffer;
        }

        public static T[] ScatterArrays<T>(Intracommunicator comm, IReadOnlyList<T[]> arraysToScatter_root, 
            int arrayLength_process, int root)
        {
            T[] allValues = null;
            int[] arrayLengths = new int[comm.Size];
            if (comm.Rank == root)
            {
                // Find out how many values there are.
                Debug.Assert(arraysToScatter_root.Count == comm.Size, 
                    $"There are {arraysToScatter_root.Count} arrays, but {comm.Size} processes."); //TODO: this will cause a deadlock in other processes
                int numValuesTotal = 0;
                foreach (T[] arry in arraysToScatter_root) numValuesTotal += arry.Length;

                // Gather all values in the same array.
                allValues = new T[numValuesTotal];
                int offset = 0;
                for (int p = 0; p < arraysToScatter_root.Count; ++p)
                {
                    T[] arry = arraysToScatter_root[p];
                    if (p == root) offset += arry.Length; // Values stored in root do no need to be processed
                    else
                    {
                        Array.Copy(arry, 0, allValues, offset, arry.Length);
                        arrayLengths[p] = arry.Length;
                        offset += arry.Length;
                    }
                }
            }
            else
            {
                // WARNING: This is a hack that depends on the weird implementation of MPI.NET needing the int[] counts array in 
                // all processes (https://github.com/mpidotnet/MPI.NET/blob/master/MPI/Intracommunicator.cs#L2008). Note that 
                // each process only accesses counts[rank], which is reasonable, since the underlying MPI implementation only 
                // needs that. 
                arrayLengths[comm.Rank] = arrayLength_process;
            }

            // Scatter the arrays
            T[] processArray = comm.ScatterFromFlattened<T>(allValues, arrayLengths, root);
            if (comm.Rank == root) processArray = arraysToScatter_root[root]; // Return the stored values instead of the scattered zeros
            return processArray;
        }

        /// <summary>
        /// Must be paired with <see cref="ReceiveArray(Intracommunicator, int, int)"/>.
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="vals"></param>
        /// <param name="dest"></param>
        /// <param name="tag"></param>
        public static void SendArray<T>(Intracommunicator comm, T[] vals, int dest, int tag)
        {
            // Sending the length is needed to allocate a buffer in destination before sending the whole array. Furthermore it 
            // blocks all processes except for the one currently processed master, preventing them from sending the whole 
            // arrays. Not sure if this is good or bad.
            comm.Send<int>(vals.Length, dest, tag);
            comm.Send<T>(vals, dest, tag);
        }
    }
}
