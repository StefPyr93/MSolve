using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;

namespace ISAAR.MSolve.Solvers.Logging
{
    public static class FetiDPDofSeparationLogging
    {
        public static void PrintDofSeparationMpi(ProcessDistribution procs, IModel model, IFetiDPDofSeparator dofSeparator)
        {
            if (procs.IsMasterProcess) Console.WriteLine("\nDof separation:");
            MpiUtilities.DoInTurn(procs.Communicator, () =>
            {
                foreach (int s in procs.GetSubdomainIDsOfProcess(procs.OwnRank))
                {
                    ISubdomain subdomain = model.GetSubdomain(s);
                    Console.Write($"Process {procs.OwnRank} - ");
                    PrintSubdomainDofSeparation(subdomain, dofSeparator);
                }
            });
        }

        public static void PrintDofSeparationSerial(IModel model, IFetiDPDofSeparator dofSeparator)
        {
            Console.WriteLine("\nDof separation:");
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                PrintSubdomainDofSeparation(subdomain, dofSeparator);
            }
        }

        private static void PrintSubdomainDofSeparation(ISubdomain subdomain, IFetiDPDofSeparator dofSeparator)
        {
            Console.WriteLine($"Subdomain {subdomain.ID}");

            // Corner dofs
            int[] cornerDofs = dofSeparator.GetCornerDofIndices(subdomain);
            Console.Write($"Corner dof indices into free dofs ({cornerDofs.Length}): ");
            WriteArray(cornerDofs);

            // Remainder dofs
            int[] remainderDofs = dofSeparator.GetRemainderDofIndices(subdomain);
            Console.Write($"Remainder dof indices into free dofs ({remainderDofs.Length}): ");
            WriteArray(remainderDofs);

            // Boundary dofs
            int[] boundaryDofs = dofSeparator.GetBoundaryDofIndices(subdomain);
            Console.Write($"Boundary dof indices into remainder dofs ({boundaryDofs.Length}): ");
            WriteArray(boundaryDofs);

            // Internal dofs
            int[] internalDofs = dofSeparator.GetInternalDofIndices(subdomain);
            Console.Write($"Internal dof indices into remainder dofs ({internalDofs.Length}): ");
            WriteArray(internalDofs);

            // Boolean mapping matrix Bc
            UnsignedBooleanMatrix Bc = dofSeparator.GetCornerBooleanMatrix(subdomain);
            Console.Write($"Bc ({Bc.NumRows} x {Bc.NumColumns}): ");
            new FullMatrixWriter() { NumericFormat = new GeneralNumericFormat() }.WriteToConsole(Bc);
        }

        private static void WriteArray<T>(T[] array)
        {
            for (int i = 0; i < array.Length; ++i) Console.Write(array[i] + " ");
            Console.WriteLine();
        }
    }
}
