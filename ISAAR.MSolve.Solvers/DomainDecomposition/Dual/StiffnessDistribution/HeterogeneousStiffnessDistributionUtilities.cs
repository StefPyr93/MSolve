using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Table<INode, IDofType, BoundaryDofLumpedStiffness> should be a dedicaed class
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    internal static class HeterogeneousStiffnessDistributionUtilities
    {
        private const int invalidDof = -1; 

        internal static double[] CalcBoundaryDofCoefficients(IFetiDPDofSeparator dofSeparator, ISubdomain subdomain, 
            Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            //TODO: Should this be cached? It stores the same info as HeterogeneousStiffnessDistribution.BoundaryDofStiffnesses.
            //      This format is more compact and has more efficient indexing when dealing with only 1 subdomain at a time, 
            //      but is difficult to index when accessing global boundary dofs. Is it possible to only use one of the two 
            //      formats? 
            //TODO: Could this be handled when extracting the lumped boundary stiffnesses? That way we can avoid searching
            //      the subdomain in BoundaryDofLumpedStiffness.SubdomainStiffnesses for each dof.

            int[] boundaryDofIndices = dofSeparator.GetBoundaryDofIndices(subdomain);
            (INode, IDofType)[] boundaryDofs = dofSeparator.GetBoundaryDofs(subdomain);
            int numBoundaryDofs = boundaryDofIndices.Length;
            var relativeStiffnesses = new double[numBoundaryDofs];
            for (int i = 0; i < boundaryDofIndices.Length; ++i)
            {
                (INode node, IDofType dofType) = boundaryDofs[i];
                relativeStiffnesses[i] = boundaryDofStiffnesses[node, dofType].CalcRelativeStiffness(subdomain);
            }
            return relativeStiffnesses;
        }

        internal static Table<INode, IDofType, BoundaryDofLumpedStiffness> CalcBoundaryDofStiffnesses(
            IFetiDPDofSeparator dofSeparator, Dictionary<ISubdomain, IIndexable2D> stiffnessesFreeFree)
        {
            var boundaryDofStiffnesses = new Table<INode, IDofType, BoundaryDofLumpedStiffness>();
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                foreach (IDofType dofType in nodeDofsPair.Value)
                {
                    var subdomainStiffnesses = new Dictionary<ISubdomain, double>();
                    double totalStiffness = 0.0;
                    foreach (ISubdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        int dofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        IIndexable2D Kff = stiffnessesFreeFree[subdomain];
                        double stiffness = Kff[dofIdx, dofIdx]; //TODO: optimized GetDiagonal(i) or GetDiagonal(subset) method for matrices. 
                        subdomainStiffnesses[subdomain] = stiffness;
                        totalStiffness += stiffness;
                    }
                    boundaryDofStiffnesses[node, dofType] = new BoundaryDofLumpedStiffness(subdomainStiffnesses, totalStiffness);
                }
            }
            return boundaryDofStiffnesses;
        }

        internal static Table<INode, IDofType, BoundaryDofLumpedStiffness> CalcBoundaryDofStiffnesses(
            IFetiDPDofSeparator dofSeparator, Dictionary<ISubdomain, IMatrixView> stiffnessesSb)
        {
            Dictionary<ISubdomain, int[]> freeToBoundaryDofmaps = MapFreeToBoundaryDofs(dofSeparator, stiffnessesSb.Keys);

            var boundaryDofStiffnesses = new Table<INode, IDofType, BoundaryDofLumpedStiffness>();
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                foreach (IDofType dofType in nodeDofsPair.Value)
                {
                    var subdomainStiffnesses = new Dictionary<ISubdomain, double>();
                    double totalStiffness = 0.0;
                    foreach (ISubdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        int freeDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        int boundaryDofIdx = freeToBoundaryDofmaps[subdomain][freeDofIdx];
                        Debug.Assert(boundaryDofIdx != invalidDof);

                        IIndexable2D Sb = stiffnessesSb[subdomain];
                        double stiffness = Sb[boundaryDofIdx, boundaryDofIdx]; //TODO: optimized GetDiagonal(i) or GetDiagonal(subset) method for matrices. 
                        subdomainStiffnesses[subdomain] = stiffness;
                        totalStiffness += stiffness;
                    }
                    boundaryDofStiffnesses[node, dofType] = new BoundaryDofLumpedStiffness(subdomainStiffnesses, totalStiffness);
                }
            }
            return boundaryDofStiffnesses;
        }

        internal static DiagonalMatrix BuildDlambda(ILagrangeMultipliersEnumerator lagrangeEnumerator,
            Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var Dlambda = new double[numLagranges];
            for (int i = 0; i < numLagranges; ++i)
            {
                LagrangeMultiplier lagrange = lagrangeEnumerator.LagrangeMultipliers[i];
                BoundaryDofLumpedStiffness boundaryDofStiffness = boundaryDofStiffnesses[lagrange.Node, lagrange.DofType];
                Dictionary<ISubdomain, double> stiffnessPerSubdomain = boundaryDofStiffness.SubdomainStiffnesses;
                double totalStiffness = boundaryDofStiffness.TotalStiffness;
                Dlambda[i] = stiffnessPerSubdomain[lagrange.SubdomainPlus] * stiffnessPerSubdomain[lagrange.SubdomainMinus]
                    / totalStiffness;
            }
            return DiagonalMatrix.CreateFromArray(Dlambda, false);
        }

        //TODO: this is also done when distributing the nodal loads. Do it here only and use the inv(Db) matrix there.
        //      Even better that code should be incorporated here, and inv(Db) should be created once and stored.
        //TODO: Kbb is also calculated for most preconditioners. Just take its diagonal and invert.
        private static DiagonalMatrix InvertBoundaryDofStiffnesses(IFetiDPDofSeparator dofSeparator, 
            ISubdomain subdomain, Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            //Debug.WriteLine($"{this.GetType().Name}: Calculating inv(diag(Kbb)) of subomain {subdomain.ID}");
            (INode node, IDofType dofType)[] boundaryDofs = dofSeparator.GetBoundaryDofs(subdomain);
            var invDb = new double[boundaryDofs.Length];
            for (int i = 0; i < boundaryDofs.Length; ++i)
            {
                (INode node, IDofType dofType) = boundaryDofs[i];
                double subdomainStiffness = boundaryDofStiffnesses[node, dofType].SubdomainStiffnesses[subdomain];
                invDb[i] = 1.0 / subdomainStiffness;
            }
            var result = DiagonalMatrix.CreateFromArray(invDb, false);
            return result;
        }

        //TODO: These should be cached somehow. Better do the whole process in dofSeparator
        private static Dictionary<ISubdomain, int[]> MapFreeToBoundaryDofs(IFetiDPDofSeparator dofSeparator, 
            IEnumerable<ISubdomain> subdomains)
        {
            var freeToBoundaryDofmaps = new Dictionary<ISubdomain, int[]>();
            {
                foreach (ISubdomain subdomain in subdomains)
                {
                    int numFreeDofs = subdomain.FreeDofOrdering.NumFreeDofs;
                    var mapFreeToBoundary = new int[numFreeDofs];
                    for (int i = 0; i < numFreeDofs; ++i) mapFreeToBoundary[i] = invalidDof;

                    int[] remainderDofs = dofSeparator.GetRemainderDofIndices(subdomain);
                    int[] boundaryDofs = dofSeparator.GetBoundaryDofIndices(subdomain);
                    for (int i = 0; i < boundaryDofs.Length; ++i)
                    {
                        int remainderIdx = boundaryDofs[i];
                        int freeIdx = remainderDofs[remainderIdx];
                        mapFreeToBoundary[freeIdx] = i;
                    }
                    freeToBoundaryDofmaps[subdomain] = mapFreeToBoundary;
                }
            }
            return freeToBoundaryDofmaps;
        }

        //TODO: This should be modified to CSR or CSC format and then benchmarked against the implicit alternative.
        /// <summary>
        /// Calculates the product Bpb = Dλ * Bb * inv(Db) explicitly, stores it and uses it for multiplications.
        /// </summary>
        public class ScalingBooleanMatrixExplicit : IMappingMatrix
        {
            private readonly Matrix explicitBpb;

            public ScalingBooleanMatrixExplicit(IFetiDPDofSeparator dofSeparator, ISubdomain subdomain, 
                Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses,
                DiagonalMatrix Dlambda, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
            {
                // According to Fragakis PhD (e.q. 3.28): 
                // Bpb = Dλ * Bb * inv(Db(s)), Dλ[λ,λ] = K(i)[b,b] * K(j)[b,b] / Sum(K(1)[b,b] + K(2)[b,b] + ...)
                // where K(s)[b,b] is the diagonal entry of (s) subdomain's stiffess matrix corresponding to the boundary dof b 
                // and (i, j) are the subdomains connected via the Lagrange multiplier λ. 
                Matrix invDb = InvertBoundaryDofStiffnesses(dofSeparator, subdomain, boundaryDofStiffnesses).CopyToFullMatrix(); //TODO: This can be reused from previous analysis steps
                this.explicitBpb = boundarySignedBooleanMatrix.MultiplyRight(invDb).MultiplyLeft(Dlambda.CopyToFullMatrix());
            }

            public int NumColumns => explicitBpb.NumColumns;

            public int NumRows => explicitBpb.NumRows;
            
            public Matrix CopyToFullMatrix() => LinearAlgebra.Commons.DenseStrategies.CopyToFullMatrix(this);

            public Vector Multiply(Vector vector, bool transposeThis = false)
                => explicitBpb.Multiply(vector, transposeThis);

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
                => explicitBpb.MultiplyRight(other, transposeThis);
        }

        /// <summary>
        /// Stores the matrices Bb, Dλ and inv(Db). Matrix-vector and matrix-matrix multiplications with Bpb = Dλ * Bb * inv(Db) 
        /// are performed implicitly, e.g. Bpb * x = Dλ * (Bb * (inv(Db) * x)).
        /// </summary>
        public class ScalingBooleanMatrixImplicit : IMappingMatrix
        {
            /// <summary>
            /// Signed boolean matrix with only the boundary dofs of the subdomain as columns. 
            /// </summary>
            private readonly SignedBooleanMatrixColMajor Bb;

            /// <summary>
            /// Diagonal matrix that stores for each dof the product of the stiffnesses corresponding to that dof in each 
            /// subdomain, divided by their sum.
            /// </summary>
            private readonly DiagonalMatrix Dlambda;

            /// <summary>
            /// Inverse of the diagonal matrix that stores the multiplicity of each boundary dof of the subdomain.
            /// </summary>
            private readonly DiagonalMatrix invDb;

            public ScalingBooleanMatrixImplicit(IFetiDPDofSeparator dofSeparator, ISubdomain subdomain,
                Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses,
                DiagonalMatrix Dlambda, SignedBooleanMatrixColMajor boundarySignedBooleanMatrix)
            {
                // According to Fragakis PhD (e.q. 3.28): 
                // Bpb = Dλ * Bb * inv(Db(s)), Dλ[λ,λ] = K(i)[b,b] * K(j)[b,b] / Sum(K(1)[b,b] + K(2)[b,b] + ...)
                // where K(s)[b,b] is the diagonal entry of (s) subdomain's stiffess matrix corresponding to the boundary dof b 
                // and (i, j) are the subdomains connected via the Lagrange multiplier λ. 
                this.Dlambda = Dlambda;
                this.Bb = boundarySignedBooleanMatrix;
                this.invDb = InvertBoundaryDofStiffnesses(dofSeparator, subdomain, boundaryDofStiffnesses); //TODO: This can be reused from previous analysis steps
            }

            public int NumColumns => invDb.NumColumns;

            public int NumRows => Dlambda.NumRows;

            public Matrix CopyToFullMatrix() => LinearAlgebra.Commons.DenseStrategies.CopyToFullMatrix(this);

            public Vector Multiply(Vector vector, bool transposeThis = false)
            {
                //TODO: Perhaps I can reuse the temporary vectors to reduce allocations/deallocations.
                if (transposeThis)
                {
                    // Bpb^T * x = (Dλ * Bb * inv(Db))^T * x = inv(Db)^T * Bb^T * Dλ^T * x = inv(Db) * (Bb^T * (Dλ * x));
                    return invDb * Bb.Multiply(Dlambda * vector, true);
                }
                else
                {
                    // Bpb * x = Dλ * (Bb * (inv(Db) * x))
                    return Dlambda * Bb.Multiply(invDb * vector);
                }
            }

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
            {
                //TODO: Perhaps I can reuse the temporary matrices to reduce allocations/deallocations.
                if (transposeThis)
                {
                    // Bpb^T * X = (Dλ * Bb * inv(Db))^T * X = inv(Db)^T * Bb^T * Dλ^T * X = inv(Db) * (Bb^T * (Dλ * X));
                    return invDb * Bb.MultiplyRight(Dlambda * other, true);
                }
                else
                {
                    // Bpb * X = Dλ * (Bb * (inv(Db) * X))
                    return Dlambda * Bb.MultiplyRight(invDb * other);
                }
            }
        }
    }
}