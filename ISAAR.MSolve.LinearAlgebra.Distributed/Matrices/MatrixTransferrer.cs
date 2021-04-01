using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using MPI;

//TODO: Repeat many of the methods of ISubdomainDataTransferrer here, but only for Matrix. These will call some of the methods
//      in ISubdomainDataTransferrer.
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Vectors
{
    public class MatrixTransferrer
    {
        private readonly ProcessDistribution procs;

        public MatrixTransferrer(ProcessDistribution processDistribution)
        {
            this.procs = processDistribution;
        }

        /// <summary>
        /// If the dimensions of the matrix are already known in each process, use 
        /// <see cref="BroadcastMatrix(ref Matrix, int, int)"/> instead
        /// instead.
        /// </summary>
        public void BroadcastMatrix(ref Matrix matrix)
        {
            var dimensions = new int[2];
            if (procs.IsMasterProcess)
            {
                dimensions[0] = matrix.NumRows;
                dimensions[1] = matrix.NumColumns;
            }
            procs.Communicator.Broadcast<int>(ref dimensions, procs.MasterProcess);
            BroadcastMatrix(ref matrix, dimensions[0], dimensions[1]);
        }

        /// <summary>
        /// More efficient than <see cref="BroadcastMatrix(ref Matrix)"/>, but each process must know the dimensions of 
        /// the matrix.
        /// </summary>
        public void BroadcastMatrix(ref Matrix matrix, int numRows, int numCols)
        {
            double[] data;
            if (procs.IsMasterProcess) data = matrix.RawData;
            else data = new double[numRows * numCols];
            procs.Communicator.Broadcast<double>(ref data, procs.MasterProcess);
            if (!procs.IsMasterProcess) matrix = Matrix.CreateFromArray(data, numRows, numCols, false);
        }

        public void SumMatrices(IEnumerable<Matrix> processMatrices, Matrix totalMatrix_master)
        {
            Matrix processSum = null;
            using (IEnumerator<Matrix> iterator = processMatrices.GetEnumerator())
            {
                bool isNotEmpty = iterator.MoveNext();
                if (!isNotEmpty) throw new ArgumentException("There must be at least 1 matrix");
                processSum = iterator.Current.Copy(); // Optimization for the first one
                while (iterator.MoveNext()) processSum.AddIntoThis(iterator.Current);
            }

            // MPI-reduce the vectors of all processes
            double[] asArray = processSum.RawData;
            double[] totalVectorData_master = null;
            if (procs.IsMasterProcess) totalVectorData_master = totalMatrix_master.RawData;
            procs.Communicator.Reduce<double>(asArray, Operation<double>.Add, procs.MasterProcess, ref totalVectorData_master);
        }
    }
}
