//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using ISAAR.MSolve.Discretization.Commons;
//using ISAAR.MSolve.Discretization.FreedomDegrees;
//using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
//using ISAAR.MSolve.LinearAlgebra.Vectors;
//using ISAAR.MSolve.Solvers.DomainDecomposition.DofSeparation;
//using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
//using ISAAR.MSolve.Solvers.LinearSystems;

//namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
//{
//    public class HeterogeneousStiffnessDistributionMpi
//    {
//        //TODO: (INode node, IDofType dofType)[] BoundaryDofs of a subdomain is stored in FetiDPDofSeparatorMpi of the 
//        //      corresponding process. HeterogeneousStiffnessDistributionMpi should gather them to master to use them.
//        //      In HomogeneousStiffnessDistributionMpi it is faster to gather the int[] arrays that result from using
//        //      (INode node, IDofType dofType)[]. Therefore they should not be always stored in master' FetiDPDofSeparatorMpi.

//        //TODO: Stiffness matrices should be accessed by LinearSystems, MatrixManagers or a dedicated class that accesses the 
//        //      correct matrices, not injected into Update()

//        //TODO: The Dλ matrix used for the preconditioner is common in all subdomains. In MPI it should be calculated by each 
//        //      process separately. In serial implementation it should only be calculated once.
//    }
//}
