using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Transfer
{
    public class DofTableSerializer
    {
        private readonly IDofSerializer dofSerializer;

        public DofTableSerializer(IDofSerializer dofSerializer)
        {
            this.dofSerializer = dofSerializer;
        }

        public int CountEntriesOf(int[] serializedTable)
        {
            Debug.Assert(serializedTable.Length % 3 == 0,
                "The provided flattened table does not correspond to a dof table. It must have a length divisible by 3, where"
                + " the first element of each triad is the node, the second is the dof type and the thrd is the dof index");
            return serializedTable.Length / 3;
        }

        public DofTable Deserialize(int[] serializedTable, Func<int, INode> getGlobalNode)
        {
            Debug.Assert(serializedTable.Length % 3 == 0,
                "The provided flattened table does not correspond to a dof table. It must have a length divisible by 3, where"
                + " the first element of each triad is the node, the second is the dof type and the thrd is the dof index");
            int numEntries = serializedTable.Length / 3;
            var table = new DofTable();
            for (int i = 0; i < numEntries; ++i)
            {
                INode node = getGlobalNode(serializedTable[3 * i]);
                IDofType dofType = dofSerializer.Deserialize(serializedTable[3 * i + 1]);
                int dofIdx = serializedTable[3 * i + 2];
                table[node, dofType] = dofIdx; //TODO: perhaps this can be optimized to avoid checking if there is already such an entry.
            }
            return table;
        }

        public int CalcPackedLength(DofTable table) => 3 * table.EntryCount;

        public void PackTableIntoArray(DofTable table, int[] buffer, int offset)
        {
            foreach ((INode row, IDofType col, int val) in table)
            {
                buffer[offset++] = row.ID;                          // Node
                buffer[offset++] = dofSerializer.Serialize(col);    // Dof type
                buffer[offset++] = val;                             // Dof index in vectors/matrices
            }
        }

        public int[] Serialize(DofTable table)
        {
            int numEntries = table.EntryCount;
            var serializedTable = new int[3 * numEntries];
            int counter = -1;
            foreach ((INode row, IDofType col, int val) in table)
            {
                serializedTable[++counter] = row.ID;                          // Node
                serializedTable[++counter] = dofSerializer.Serialize(col);    // Dof type
                serializedTable[++counter] = val;                             // Dof index in vectors/matrices
            }
            return serializedTable;
        }

        public DofTable UnpackTableFromArray(int[] buffer, int start, int end, Func<int, INode> getGlobalNode)
        {
            int numEntries = (end - start) / 3;
            var table = new DofTable();
            for (int i = 0; i < numEntries; ++i)
            {
                INode node = getGlobalNode(buffer[start + 3 * i]);
                IDofType dofType = dofSerializer.Deserialize(buffer[start + 3 * i + 1]);
                int dofIdx = buffer[start + 3 * i + 2];
                table[node, dofType] = dofIdx; //TODO: perhaps this can be optimized to avoid checking if there is already such an entry.
            }
            return table;
        }
    }
}
