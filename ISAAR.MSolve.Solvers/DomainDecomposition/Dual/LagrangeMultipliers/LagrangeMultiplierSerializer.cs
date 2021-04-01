using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;

//TODO: The dictionary of subdomain nodes should be accessed by the subdomain.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    public class LagrangeMultiplierSerializer
    {
        private readonly IDofSerializer dofSerializer;

        public LagrangeMultiplierSerializer(IDofSerializer dofSerializer)
        {
            this.dofSerializer = dofSerializer;
        }

        public (int numGlobalLagranges, List<SubdomainLagrangeMultiplier> subdomainLagranges) Deserialize(
            int[] serializedLagranges, ISubdomain subdomain)
        {
            CheckSerializedLength(serializedLagranges);
            int numGlobalLagranges = serializedLagranges.Length / 4;
            var subdomainLagranges = new List<SubdomainLagrangeMultiplier>();
            for (int i = 0; i < numGlobalLagranges; ++i)
            {
                if (serializedLagranges[4 * i + 2] == subdomain.ID)
                {
                    INode node = subdomain.GetNode(serializedLagranges[4 * i]);
                    IDofType dofType = dofSerializer.Deserialize(serializedLagranges[4 * i + 1]);
                    subdomainLagranges.Add(new SubdomainLagrangeMultiplier(i, node, dofType, true));
                }
                else if (serializedLagranges[4 * i + 3] == subdomain.ID)
                {
                    INode node = subdomain.GetNode(serializedLagranges[4 * i]);
                    IDofType dofType = dofSerializer.Deserialize(serializedLagranges[4 * i + 1]);
                    subdomainLagranges.Add(new SubdomainLagrangeMultiplier(i, node, dofType, false));
                }
            }
            return (numGlobalLagranges, subdomainLagranges);
        }

        //TODO: Not thrilled about having an array of DTOs with null entries in the array and the DTOs. Especially if it is passed between classes
        /// <summary>
        /// For the lagrange multipliers that are applied to the subdomain corresponding to this process, the opposite subdomain
        /// may be null. For all other lagrange multipliers, the corresponding array entries may be null. 
        /// </summary>
        /// <param name="serializedLagranges"></param>
        /// <param name="subdomain"></param>
        public LagrangeMultiplier[] DeserializeIncompletely(int[] serializedLagranges,
            Dictionary<int, ISubdomain> processSubdomains)
        {
            CheckSerializedLength(serializedLagranges);
            int numLagranges = serializedLagranges.Length / 4;
            var lagranges = new LagrangeMultiplier[numLagranges];
            for (int i = 0; i < numLagranges; ++i)
            {
                int subdomainPlusID = serializedLagranges[4 * i + 2];
                int subdomainMinusID = serializedLagranges[4 * i + 3];
                bool subdomainPlusIsStored = processSubdomains.TryGetValue(subdomainPlusID, out ISubdomain subdomainPlus);
                bool subdomainMinusIsStored = processSubdomains.TryGetValue(subdomainMinusID, out ISubdomain subdomainMinus);

                if (subdomainPlusIsStored || subdomainMinusIsStored)
                {
                    int nodeID = serializedLagranges[4 * i];
                    INode node = subdomainPlusIsStored ? subdomainPlus.GetNode(nodeID) : subdomainMinus.GetNode(nodeID);
                    IDofType dofType = dofSerializer.Deserialize(serializedLagranges[4 * i + 1]);
                    lagranges[i] = new LagrangeMultiplier(node, dofType, subdomainPlus, subdomainMinus); // one of them may be null
                }
            }
            return lagranges;
        }

        public int[] Serialize(IReadOnlyList<LagrangeMultiplier> lagranges)
        {
            var serializedLagranges = new int[4 * lagranges.Count];
            for (int i = 0; i < lagranges.Count; ++i)
            {
                LagrangeMultiplier lagr = lagranges[i];
                serializedLagranges[4 * i] = lagr.Node.ID;
                serializedLagranges[4 * i + 1] = dofSerializer.Serialize(lagr.DofType);
                serializedLagranges[4 * i + 2] = lagr.SubdomainPlus.ID;
                serializedLagranges[4 * i + 3] = lagr.SubdomainMinus.ID;
            }
            return serializedLagranges;
        }

        [Conditional("DEBUG")]
        private void CheckSerializedLength(int[] serializedLagranges)
        {
            if (serializedLagranges.Length % 4 != 0)
            {
                throw new ArgumentException("The provided int[] array of serialized lagrange multipliers is not valid."
                + " The array's length must be divisible by 4: element 0 = node, element 1 = dof type,"
                + " element 2 = plus subdomain, element 3 = minus subdomain.");
            }
        }
    }
}
