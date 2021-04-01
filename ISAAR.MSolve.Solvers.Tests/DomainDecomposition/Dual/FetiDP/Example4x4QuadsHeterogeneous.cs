using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;

//TODO: Split this into Homogeneous and Heterogeneous classes. Everything is different but the connectivity.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Example4x4QuadsHeterogeneous
    {
        public static Vector SolutionGlobalDisplacements => Vector.CreateFromArray(new double[]
        {
            17.623494584618864, 12.564560593215612, 31.832863897135404, 34.496634608059082, 40.255481382985629,
            66.49190654178912, 42.572002358887204, 99.798764204232072, 4.267568672307144, 9.00506902466324,
            9.100928263505315, 31.107370029452451, 12.1615036308774, 66.065492717632239, 11.510673148931499,
            102.06649895017948, -3.0529124682202156, 9.24107474483673, -7.8531777412741217, 26.728892403726846,
            -16.890006178831449, 70.602493468916791, -21.80233265288679, 109.39882637058051, -4.7311061272016808,
            10.030926199331375, -5.6722429958962142, 18.837815470700932, 146.94209278892487, 392.04674590737193,
            -35.619167413693908, 1407.200332011206, -9.9609496807814057, 10.46574373452243, -17.603838651152756,
            20.760800663270086, -843.13592713307355, 371.10700308359418, -1666.2547486301742, 3714.1637893447919
        });

        public static Model CreateModel(double stiffnessRatio)
        {
            //                                    Λ P
            //                                    | 
            //                                     
            // |> 20 ---- 21 ---- 22 ---- 23 ---- 24
            //    |  (12) |  (13) |  (14) |  (15) |
            //    |  E0   |  E0   |  E1   |  E1   |
            // |> 15 ---- 16 ---- 17 ---- 18 ---- 19
            //    |  (8)  |  (9)  |  (10) |  (11) |
            //    |  E0   |  E0   |  E1   |  E1   |
            // |> 10 ---- 11 ---- 12 ---- 13 ---- 14
            //    |  (4)  |  (5)  |  (6)  |  (7)  |
            //    |  E0   |  E0   |  E0   |  E0   |
            // |> 5 ----- 6 ----- 7 ----- 8 ----- 9
            //    |  (0)  |  (1)  |  (2)  |  (3)  |
            //    |  E0   |  E0   |  E0   |  E0   |
            // |> 0 ----- 1 ----- 2 ----- 3 ----- 4


            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 4.0;
            builder.DomainLengthY = 4.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 4;
            builder.NumTotalElementsY = 4;
            //builder.YoungModulus = 1.0;
            double E0 = 1.0;
            builder.YoungModuliOfSubdomains = new double[,]
            {
                { E0, E0 }, { E0, stiffnessRatio * E0}
            };
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.UpperRightCorner, StructuralDof.TranslationY, 10.0);

            return builder.BuildModel();
        }

        public static UsedDefinedCornerNodes DefineCornerNodeSelectionSerial(Model model) 
            => Example4x4QuadsHomogeneous.DefineCornerNodeSelectionSerial(model);

        public static Matrix GetMatrixBpbr(int subdomainID)
        {
            if (subdomainID == 0)
            {
                return Matrix.CreateFromArray(new double[,]
                {
                    { 0.5, 0, 0, 0 },
                    { 0, 0.5, 0, 0 },
                    { 0, 0, 0.5, 0 },
                    { 0, 0, 0, 0.5 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 }
                });
            }
            else if (subdomainID == 1)
            {
                return Matrix.CreateFromArray(new double[,]
                {
                    { -0.5, 0, 0, 0 },
                    { 0, -0.5, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0.00990099009900990, 0 },
                    { 0, 0, 0, 0.00990099009900990 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 }
                });
            }
            else if (subdomainID == 2)
            {
                return Matrix.CreateFromArray(new double[,]
                {
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { -0.5, 0, 0, 0 },
                    { 0, -0.5, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0.00990099009900990, 0 },
                    { 0, 0, 0, 0.00990099009900990 }
                });
            }
            else if (subdomainID == 3)
            {
                return Matrix.CreateFromArray(new double[,]
                {
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { 0, 0, 0, 0 },
                    { -0.990099009900990, 0, 0, 0 },
                    { 0, -0.990099009900990, 0, 0 },
                    { 0, 0, -0.990099009900990, 0 },
                    { 0, 0, 0, -0.990099009900990 }
                });
            }
            else throw new ArgumentException("Subdomain ID must be 0, 1, 2 or 3");
        }
    }
}
