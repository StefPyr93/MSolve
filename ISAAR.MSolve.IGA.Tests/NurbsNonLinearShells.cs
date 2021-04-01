using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
//using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
//using ISAAR.MSolve.MultiscaleAnalysis;
//using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Troschuetz.Random;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsNonLinearShells
    {
        private const double Tolerance = 1e-9;

        #region Fields

        private int degreeKsi = 2;
        private int degreeHeta = 2;
        private List<ControlPoint> ElementControlPoints()
        {
            return new List<ControlPoint>
            {
                new ControlPoint {ID = 0, X = 0.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 1, X = 0.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 2, X = 0.0, Y =  1, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 3, X = 0.3125, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 4, X = 0.3125, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 5, X = 0.3125, Y =  1, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 6, X = 0.9375, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 7, X = 0.9375, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 8, X = 0.9375, Y =  1, Z = 0.0, WeightFactor = 1.0},
            };
        }

        private List<Knot> ElementKnots()
        {
            return new List<Knot>()
            {
                new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=1,Ksi=0.0,Heta=1.0,Zeta =0.0 },
                new Knot(){ID=2,Ksi=0.625,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=3,Ksi=0.625,Heta=1,Zeta =0.0 }
            };
        }

        private double[] knotValueVectorKsi=>new double[]
        {
            0.0000000000000, 0.0000000000000, 0.0000000000000, 0.6250000000000, 1.2500000000000, 1.8750000000000,
            2.5000000000000, 3.1250000000000, 3.7500000000000, 4.3750000000000, 5.0000000000000, 5.6250000000000,
            6.2500000000000, 6.8750000000000, 7.5000000000000, 8.1250000000000, 8.7500000000000, 9.3750000000000,
            10.0000000000000, 10.0000000000000, 10.0000000000000,
        };
        

        private double[] knotValueVectorHeta=>new double[]
        {
            0.0, 0.0, 0.0, 1.0, 1.0, 1.0
        };
    
        private ShellElasticMaterial2D Material => new ShellElasticMaterial2D()
        {
            YoungModulus = 1200000,
            PoissonRatio = 0.0,
            TangentVectorV1 = new double[] { 1.000000000000000000000000000000000000, 0.000000000000000053342746886286800000, 0.000000000000000000000000000000000000 },
            TangentVectorV2 = new double[] { 3.90312782094781000000000000000000E-18, 9.99999999999999000000000000000000E-01, 0.00000000000000000000000000000000E+00, },
            NormalVectorV3 = new double[] { 0, 0, 1 }
        };
        
        private double[] localSolution => new double[]
        {
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            -0.000001762030370944800000000000000000000, 0.000000000000000018440454162078200000000,
            0.001664779366416430000000000000000000000, -0.000001762030370985190000000000000000000,
            0.000000000000000001055844403292750000000, 0.001664779366416240000000000000000000000,
            -0.000001762030370843590000000000000000000, -0.000000000000000079088297192718400000000,
            0.001664779366416050000000000000000000000
        };

        private Forces MembraneForces => new Forces()
        {
            v0 = -0.0327209647710278000000000000000000000000000000000,
            v1 = -0.0000000000000000000000000000193609885449360000000,
            v2 = 0.0000000000001690940892323080000000000000000000000
        };

        private Forces BendingMoments => new Forces()
        {
            v0 = -0.42618363401210900000000000000000000000000000000000,
            v1 = -0.00000000000000000075669834331405100000000000000000,
            v2 = 0.00000000000000683985998006887000000000000000000000,
        };

            #endregion

        [Fact]
        public void JacobianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = new double[2, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var expectedJacobian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "jacobian");
            
            for (var i = 0; i < jacobianMatrix.GetLength(0); i++)
            {
                for (var j = 0; j < jacobianMatrix.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedJacobian[i, j], jacobianMatrix[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void HessianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var expectedHessian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "hessian");

            for (var i = 0; i < hessian.GetLength(0); i++)
            {
                for (var j = 0; j < hessian.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedHessian[i, j], hessian[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void MembraneDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);

            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;
            var Bmembrane = new double[3, controlPoints.Length * 3];
            shellElement.CalculateMembraneDeformationMatrix(controlPoints.Length, nurbs, 0, surfaceBasisVector1, surfaceBasisVector2, Bmembrane);

            var expectedBmembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bmembrane");

            for (var i = 0; i < Bmembrane.GetLength(0); i++)
            {
                for (var j = 0; j < Bmembrane.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBmembrane[i, j], Bmembrane[i, j],
                        Tolerance));
                }
            }
        }
        
        [Fact]
        public void BendingDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);

            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector(hessian, 2);
            var Bbending = new double[3, controlPoints.Length * 3];
            shellElement.CalculateBendingDeformationMatrix(controlPoints.Length, surfaceBasisVector3, nurbs, 0, surfaceBasisVector2,
                surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, Bbending);

            var expectedBbending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bbending");

            for (var i = 0; i < Bbending.GetLength(0); i++)
            {
                for (var j = 0; j < Bbending.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBbending[i, j], Bbending[i, j],
                        Tolerance));
                }
            }
        }
        
        
        [Fact]
        public void StiffnessMatrixMembraneNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);

            var KmembraneNL = new double[controlPoints.Length * 3, controlPoints.Length * 3];

            var membraneForces = new Forces()
            {
                v0 = -0.0327209647710278000000000000000000000000000000000,
                v1 = -0.0000000000000000000000000000193609885449360000000,
                v2 = 0.0000000000001690940892323080000000000000000000000
            };

            shellElement.CalculateKmembraneNL(controlPoints, ref membraneForces, nurbs, 0, KmembraneNL);

            var expectedKmembraneNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KmembraneNL");

            for (var i = 0; i < KmembraneNL.GetLength(0); i++)
            {
                for (var j = 0; j < KmembraneNL.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKmembraneNL[i, j], KmembraneNL[i, j], Tolerance));
                }
            }
        }

        [Fact]
        public void StiffnessMatrixBendingNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                gaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);

            shellElement._solution = localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector(hessian, 2);

            var KbendingNL = new double[elementControlPoints.Length * 3, elementControlPoints.Length * 3];

            var BendingMoments = new Forces()
            {
                v0 = -0.42618363401210900000000000000000000000000000000000,
                v1 = -0.00000000000000000075669834331405100000000000000000,
                v2 = 0.00000000000000683985998006887000000000000000000000,
            };

            shellElement.CalculateKbendingNL(elementControlPoints, ref BendingMoments, nurbs,
                surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3,
                surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, J1, 0, KbendingNL);

            var expectedKbendingNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KbendingNL");


            for (var i = 0; i < KbendingNL.GetLength(0); i++)
            {
                for (var j = 0; j < KbendingNL.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKbendingNL[i, j], KbendingNL[i, j],1e-04));
                }
            }
        }
        
        [Fact]
        public void ConstitutiveMatrixThicknessIntegration()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var midsurfaceGaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                midsurfaceGaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);

            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();

            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
                shellElement.IntegratedConstitutiveOverThickness(gaussPoints[0]);

            var expectedConstitutiveMembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "MembraneConstitutiveMatrix");
            var expectedConstitutiveBending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "BendingConstitutiveMatrix");

            for (int i = 0; i < MembraneConstitutiveMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < MembraneConstitutiveMatrix.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveMembrane[i, j], MembraneConstitutiveMatrix[i, j], Tolerance));
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveBending[i, j], BendingConstitutiveMatrix[i, j], Tolerance));
                }
            }
        }

        [Fact]
        public void StressesThicknessIntegration()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var midsurfaceGaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                midsurfaceGaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            foreach (var controlPoint in controlPoints)
            {
                shellElement.ControlPointsDictionary.Add(controlPoint.ID,controlPoint);
            }

            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);

            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var MembraneForces = new Forces();
            var BendingMoments = new Forces();
            shellElement.IntegratedStressesOverThickness(gaussPoints[0], ref MembraneForces, ref BendingMoments);

            var expectedMembraneForces = new double[3]
            {
                -0.03272096477102780000000000000000000000000000000,
                -0.00000000000000000000000000001936098854493600000,
                0.00000000000016909408923230800000000000000000000,
            };

            var expectedBendingMoments = new double[3]
            {
                -0.426183634012109000000000000000000000000,
                -0.000000000000000000756698343314051000000,
                0.000000000000006839859980068870000000000,

            };

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[0], MembraneForces.v0, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[0], BendingMoments.v0, Tolerance));

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[1], MembraneForces.v1, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[1], BendingMoments.v1, Tolerance));

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[2], MembraneForces.v2, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[2], BendingMoments.v2, Tolerance));
        }

        [Fact]
        public void TotalElementStiffnessMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var midsurfaceGaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                midsurfaceGaussPoints);

            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            foreach (var controlPoint in controlPoints)
            {
                shellElement.ControlPointsDictionary.Add(controlPoint.ID,controlPoint);
            }
            
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var Ktotal=shellElement.StiffnessMatrix(shellElement);

            var expectedKtotal = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Ktotal");
            
            for (var i = 0; i < Ktotal.NumRows; i++)
            {
                for (var j = 0; j < Ktotal.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKtotal[i, j], Ktotal[i, j],
                        1e-6));
                }
            }
        }

        [Fact]
        public void InternalForcesTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var gauss= new GaussQuadrature();
            var midsurfaceGaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnots().ToArray()).ToArray();
            var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta, controlPoints,
                midsurfaceGaussPoints);
            var shellElement = new KirchhoffLoveShellNL(Material, ElementKnots(), nurbs, ElementControlPoints(),
                new Patch(), 0.1, degreeKsi, degreeHeta);
            shellElement._solution= localSolution;
            foreach (var controlPoint in controlPoints)
            {
                shellElement.ControlPointsDictionary.Add(controlPoint.ID,controlPoint);
            }

            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var f = shellElement.CalculateForces(shellElement,localSolution, new double[27]);

            var expectedForces = new double[27]
            {
                0.01333444930133850000000000000000000,
                -0.00000000000015109302777189900000000,
                0.45460748070705600000000000000000000,
                0.01333444930010990000000000000000000,
                0.00000000000029570495208513400000000,
                0.45460748070707000000000000000000000,
                0.01333444930077590000000000000000000,
                0.00000000000096541364505914100000000,
                0.45460748070708400000000000000000000,
                -0.01091686339304530000000000000000000,
                0.00000000000008994291058712290000000,
                -0.68190542227462800000000000000000000,
                -0.01091686339111670000000000000000000,
                -0.00000000000009831769021815550000000,
                -0.68190542227458800000000000000000000,
                -0.01091686339457100000000000000000000,
                -0.00000000000024966598145038100000000,
                -0.68190542227455400000000000000000000,
                -0.00241758590824831000000000000000000,
                0.00000000000006115015407641720000000,
                0.22729794156757200000000000000000000,
                -0.00241758590732423000000000000000000,
                -0.00000000000019738726144919400000000,
                0.22729794156752000000000000000000000,
                -0.00241758590791879000000000000000000,
                -0.00000000000071574770178168100000000,
                0.22729794156746600000000000000000000,
            };

            for (var j = 0; j < 27; j++)
            {
                Assert.True(Utilities.AreValuesEqual(expectedForces[j], f[j],
                    Tolerance));
            }
        }

        [Fact]
        public void TestA3r()
        {
            var surfaceBasisVector1 = new double[] { 0.49281343102090647, -0.0027907619686642761, -2.0647413184286493E-07 };
            var surfaceBasisVector2 = new double[] {0.0032374497982368571, 0.57141938776238976, 2.3218349204021492E-09};
            var surfaceBasisVector3 = new double[] { 4.1893372881671206E-07, -6.4367991618857595E-09, 0.99999999999991229};

            double dksi_r=-2.0934480533670992;
            double dheta_r=-2.0934480533670996;
            double J1=0.28161218398684618;
            var a3r= new a3r();

            KirchhoffLoveShellNL.CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r, dheta_r, J1, ref a3r);
             
            Assert.Equal(1.7882446740031873E-06,a3r.a3r00);
            Assert.Equal(-2.7475877512613432E-08,a3r.a3r01);
            Assert.Equal(4.2685621877569231,a3r.a3r02);
                         
            Assert.Equal(1.5246711353875095E-06,a3r.a3r10);
            Assert.Equal(-2.3426144068498867E-08,a3r.a3r11);
            Assert.Equal(3.639408886207,a3r.a3r12);

            Assert.Equal(-7.38802532863679E-13,a3r.a3r20);
            Assert.Equal(1.1827148338265939E-14,a3r.a3r21);
            Assert.Equal(-1.7648185299346879E-06,a3r.a3r22);
        }

        [Fact]
        public void TestCalculateCrossProduct()
        {
            var v1 = new double[] {1, 20, 50};
            var v2 = new double[] {4, 30, 45};
            var c = new double[3];
            KirchhoffLoveShellNL.CalculateCrossProduct(v1, v2, c);

            var vector1 = Vector.CreateFromArray(v1);
            var vector2 = Vector.CreateFromArray(v2);

            var cross = vector1.CrossProduct(vector2);

            Assert.Equal(cross[0],c[0]);
            Assert.Equal(cross[1],c[1]);
            Assert.Equal(cross[2],c[2]);
        }

        [Fact]
        public void TestBab_rs()
        {
            #region Initialization
            var surfaceBasisVectorDerivative1 = new double[3]
            {
                -0.0002215697315227748,
                -0.040236627065225024,
                -3.1283539134440684E-06
            };
            var surfaceBasisVectorDerivative2 = new double[3]
            {
                4.3550534925664095E-11,
                1.0323258744575857E-08,
                -1.0202863489210098E-09
            };
            var surfaceBasisVectorDerivative12 = new double[3]
            {
                0.046627259953136886,
                -0.00026442650364636821,
                6.5202798291920747E-08
            };

            var d2ksi_dr2 = 4.4992901171737891;
            var d2ksi_ds2 = 4.4992901171737891;

            var d2heta_dr2 = 4.49929011717379;
            var d2heta_ds2 = 4.49929011717379;

            var d2ksiheta_dr2 = 6.7489351757606837;
            var d2ksiheta_ds2 = 6.7489351757606837;
            var a3rs= new a3rs();
            a3rs.a3rs00_0=	1.5266467195814444E-05	;
            a3rs.a3rs00_1=	1.3016307114560258E-05	;
            a3rs.a3rs00_2=	-1.1837641977763269E-11	;
            a3rs.a3rs01_0=	6.390871065454354E-06	;
            a3rs.a3rs01_1=	5.448905725897356E-06	;
            a3rs.a3rs01_2=	-2.55440113505756E-12	;
            a3rs.a3rs02_0=	18.220623150741094	;
            a3rs.a3rs02_1=	15.535043157448497	;
            a3rs.a3rs02_2=	-2.0715372921711805E-05	;
            a3rs.a3rs10_0=	6.390871065454354E-06	;
            a3rs.a3rs10_1=	5.448905725897356E-06	;
            a3rs.a3rs10_2=	-2.55440113505756E-12	;
            a3rs.a3rs11_0=	-1.9999190555149591E-07	;
            a3rs.a3rs11_1=	-1.7051463378493538E-07	;
            a3rs.a3rs11_2=	8.8817841970012523E-14	;
            a3rs.a3rs12_0=	15.53504315745124	;
            a3rs.a3rs12_1=	13.245297041003683	;
            a3rs.a3rs12_2=	-6.22035643166942E-06	;
            a3rs.a3rs20_0=	18.220623150741094	;
            a3rs.a3rs20_1=	15.535043157448497	;
            a3rs.a3rs20_2=	-2.0715372921711805E-05	;
            a3rs.a3rs21_0=	15.53504315745124	;
            a3rs.a3rs21_1=	13.245297041003683	;
            a3rs.a3rs21_2=	-6.22035643166942E-06	;
            a3rs.a3rs22_0=	-2.8248610566845741E-05	;
            a3rs.a3rs22_1=	-1.2643252672057042E-05	;
            a3rs.a3rs22_2=	-31.465920191744779	;

            var a3r= new a3r();
            a3r.a3r00=	1.7882446740031873E-06  ;
            a3r.a3r01=	-2.7475877512613432E-08	;
            a3r.a3r02=	4.2685621877569231	;
            a3r.a3r10=	1.5246711353875095E-06	;
            a3r.a3r11=	-2.3426144068498867E-08	;
            a3r.a3r12=	3.639408886207	;
            a3r.a3r20=	-7.38802532863679E-13	;
            a3r.a3r21=	1.1827148338265939E-14	;
            a3r.a3r22=	-1.7648185299346879E-06	;

            var a3s= new a3r();
            a3s.a3r00=1.7882446740031873E-06	;
            a3s.a3r01=-2.7475877512613432E-08	;
            a3s.a3r02=4.2685621877569231	;
            a3s.a3r10=1.5246711353875095E-06	;
            a3s.a3r11=-2.3426144068498867E-08	;
            a3s.a3r12=3.639408886207	;
            a3s.a3r20=-7.38802532863679E-13	;
            a3s.a3r21=1.1827148338265939E-14	;
            a3s.a3r22=	-1.7648185299346879E-06	;
            
            var da3_drds = new double[3,3][];
            da3_drds[0, 0] = new double[]
            {
                1.5266467195814444E-05,
                1.3016307114560258E-05,
                -1.1837641977763269E-11,
            };
            da3_drds[0,1] = new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            };
            da3_drds[0, 2] = new double[]
            {
                18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05
            };

            da3_drds[1, 0] = new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            };
            da3_drds[1, 1] = new double[]
            {
                -1.9999190555149591E-07,
                -1.7051463378493538E-07,
                8.8817841970012523E-14
            };
            da3_drds[1, 2] = new double[]
            {
                15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06
            };

            da3_drds[2, 0] = new double[]
            {
                18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05
            };
            da3_drds[2, 1] = new double[]
            {
                15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06
            };
            da3_drds[2, 2] = new double[]
            {
                -2.8248610566845741E-05,
                -1.2643252672057042E-05,
                -31.465920191744779
            };

            Bab_rs Bab_rsExpected= new Bab_rs();
            Bab_rsExpected.Bab_rs00_0=	1.5564548295526568E-05	;
            Bab_rsExpected.Bab_rs00_1=	1.6091663312697994E-05	;
            Bab_rsExpected.Bab_rs00_2=	4.9691772888834864E-05	;
            Bab_rsExpected.Bab_rs01_0=	6.5156542160513032E-06	;
            Bab_rsExpected.Bab_rs01_1=	6.7363158837647762E-06	;
            Bab_rsExpected.Bab_rs01_2=	2.0802043424519865E-05	;
            Bab_rsExpected.Bab_rs02_0=	18.576384789429813	;
            Bab_rsExpected.Bab_rs02_1=	19.205499827078938	;
            Bab_rsExpected.Bab_rs02_2=	59.307438707759943	;
            Bab_rsExpected.Bab_rs10_0=	6.5156542160513032E-06	;
            Bab_rsExpected.Bab_rs10_1=	6.7363158837647762E-06	;
            Bab_rsExpected.Bab_rs10_2=	2.0802043424519861E-05	;
            Bab_rsExpected.Bab_rs11_0=	-2.0389679110046287E-07	;
            Bab_rsExpected.Bab_rs11_1=	-2.1080203875074925E-07	;
            Bab_rsExpected.Bab_rs11_2=	-6.5096608290578748E-07	;
            Bab_rsExpected.Bab_rs12_0=	15.838368261336552	;
            Bab_rsExpected.Bab_rs12_1=	16.374756571476876	;
            Bab_rsExpected.Bab_rs12_2=	50.565977478394949	;
            Bab_rsExpected.Bab_rs20_0=	18.576384789429813	;
            Bab_rsExpected.Bab_rs20_1=	19.205499827078938	;
            Bab_rsExpected.Bab_rs20_2=	59.307438707759943	;
            Bab_rsExpected.Bab_rs21_0=	15.838368261336552	;
            Bab_rsExpected.Bab_rs21_1=	16.374756571476876	;
            Bab_rsExpected.Bab_rs21_2=	50.565977478394956	;
            Bab_rsExpected.Bab_rs22_0=	8.3070654310999034E-05	;
            Bab_rsExpected.Bab_rs22_1=	-1.5848757023602572E-05	;
            Bab_rsExpected.Bab_rs22_2=	-5.4373539710938839E-05	;

            #endregion
            
            Bab_rs Bab_rsAlternative = KirchhoffLoveShellNL.CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12,d2ksi_dr2,d2ksi_ds2, d2heta_dr2,d2heta_ds2, d2ksiheta_dr2,d2ksiheta_ds2,
                a3rs, a3r, a3s, da3_drds);


            Assert.Equal(Bab_rsExpected.Bab_rs00_0,Bab_rsAlternative.Bab_rs00_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs00_1,Bab_rsAlternative.Bab_rs00_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs00_2,Bab_rsAlternative.Bab_rs00_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_0,Bab_rsAlternative.Bab_rs01_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_1,Bab_rsAlternative.Bab_rs01_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_2,Bab_rsAlternative.Bab_rs01_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_0,Bab_rsAlternative.Bab_rs02_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_1,Bab_rsAlternative.Bab_rs02_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_2,Bab_rsAlternative.Bab_rs02_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_0,Bab_rsAlternative.Bab_rs10_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_1,Bab_rsAlternative.Bab_rs10_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_2,Bab_rsAlternative.Bab_rs10_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_0,Bab_rsAlternative.Bab_rs11_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_1,Bab_rsAlternative.Bab_rs11_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_2,Bab_rsAlternative.Bab_rs11_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_0,Bab_rsAlternative.Bab_rs12_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_1,Bab_rsAlternative.Bab_rs12_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_2,Bab_rsAlternative.Bab_rs12_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_0,Bab_rsAlternative.Bab_rs20_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_1,Bab_rsAlternative.Bab_rs20_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_2,Bab_rsAlternative.Bab_rs20_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_0,Bab_rsAlternative.Bab_rs21_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_1,Bab_rsAlternative.Bab_rs21_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_2,Bab_rsAlternative.Bab_rs21_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_0,Bab_rsAlternative.Bab_rs22_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_1,Bab_rsAlternative.Bab_rs22_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_2,Bab_rsAlternative.Bab_rs22_2  ,9);
        }

        [Fact]
        public void Test_a3rs()
        {
            var surfaceBasisVector1 = Vector.CreateFromArray(new double[]{0.49281343102090647,
                                                                          -0.0027907619686642761,
                                                                          -2.0647413184286493E-07});
            var surfaceBasisVector2 =Vector.CreateFromArray(new double[]{0.0032374497982368571,
                0.57141938776238976,
                2.3218349204021492E-09});
            var surfaceBasisVector3 = Vector.CreateFromArray(new double[]{4.1893372881671206E-07,
                                                                          -6.4367991618857595E-09,
                0.99999999999991229});
            double J1 = 0.28161218398684618;
            var dksi_r = -2.0934480533670992;
            var dheta_r = -2.0934480533670996;

            var dksi_s = -2.0934480533670992;
            var dheta_s= -2.0934480533670996;

            #region Expected Values
            var a3rsExpected= new a3rs();
            a3rsExpected.a3rs00_0=	1.5266467195814444E-05	;
            a3rsExpected.a3rs00_1=	1.3016307114560258E-05	;
            a3rsExpected.a3rs00_2=	-1.1837641977763269E-11	;
            a3rsExpected.a3rs01_0=	6.390871065454354E-06	;
            a3rsExpected.a3rs01_1=	5.448905725897356E-06	;
            a3rsExpected.a3rs01_2=	-2.55440113505756E-12	;
            a3rsExpected.a3rs02_0=	18.220623150741094	;
            a3rsExpected.a3rs02_1=	15.535043157448497	;
            a3rsExpected.a3rs02_2=	-2.0715372921711805E-05	;
            a3rsExpected.a3rs10_0=	6.390871065454354E-06	;
            a3rsExpected.a3rs10_1=	5.448905725897356E-06	;
            a3rsExpected.a3rs10_2=	-2.55440113505756E-12	;
            a3rsExpected.a3rs11_0=	-1.9999190555149591E-07	;
            a3rsExpected.a3rs11_1=	-1.7051463378493538E-07	;
            a3rsExpected.a3rs11_2=	8.8817841970012523E-14	;
            a3rsExpected.a3rs12_0=	15.53504315745124	;
            a3rsExpected.a3rs12_1=	13.245297041003683	;
            a3rsExpected.a3rs12_2=	-6.22035643166942E-06	;
            a3rsExpected.a3rs20_0=	18.220623150741094	; 
            a3rsExpected.a3rs20_1=	15.535043157448497	;
            a3rsExpected.a3rs20_2=	-2.0715372921711805E-05	;
            a3rsExpected.a3rs21_0=	15.53504315745124	;
            a3rsExpected.a3rs21_1=	13.245297041003683	;
            a3rsExpected.a3rs21_2=	-6.22035643166942E-06	;
            a3rsExpected.a3rs22_0=	-2.8248610566845741E-05	;
            a3rsExpected.a3rs22_1=	-1.2643252672057042E-05	;
            a3rsExpected.a3rs22_2=	-31.465920191744779	;

            var da3_drdsExpected = new Vector[3, 3];
            da3_drdsExpected[0, 0] = Vector.CreateFromArray(new double[] { 1.5266467195814444E-05,
                1.3016307114560258E-05,
                -1.1837641977763269E-11 });
            da3_drdsExpected[0, 1] = Vector.CreateFromArray(new double[] { 6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12 });
            da3_drdsExpected[0, 2] = Vector.CreateFromArray(new double[] {18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05 });

            da3_drdsExpected[1, 0] = Vector.CreateFromArray(new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            });
            da3_drdsExpected[1, 1] = Vector.CreateFromArray(new double[] { -1.9999190555149591E-07,
                                                                           -1.7051463378493538E-07,
                8.8817841970012523E-14 });
            da3_drdsExpected[1, 2] = Vector.CreateFromArray(new double[] { 15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06 });

            da3_drdsExpected[2, 0] = Vector.CreateFromArray(new double[] { 18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05 });
            da3_drdsExpected[2, 1] = Vector.CreateFromArray(new double[] { 15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06 });
            da3_drdsExpected[2, 2] = Vector.CreateFromArray(new double[] { -2.8248610566845741E-05,
                                                                           -1.2643252672057042E-05,
                                                                           -31.465920191744779 });
            #endregion

            (a3rs a3rsAlternative, var da3tilde_drds, var da3tilde_dr, var da3tilde_ds,
                    double[] dnorma3_dr, double[] dnorma3_ds, double[,] dnorma3_drds, var a3_tilde, var da3_drds) =
                KirchhoffLoveShellNL.Calculate_a3rs(surfaceBasisVector1, surfaceBasisVector2,
                    surfaceBasisVector3, J1, dksi_r, dksi_s, dheta_r,dheta_s);
            
            #region Assertions
            Assert.Equal(a3rsExpected.a3rs00_0 ,a3rsAlternative.a3rs00_0,9);
            Assert.Equal(a3rsExpected.a3rs00_1 ,a3rsAlternative.a3rs00_1,9);
            Assert.Equal(a3rsExpected.a3rs00_2 ,a3rsAlternative.a3rs00_2,9);
            Assert.Equal(a3rsExpected.a3rs01_0 ,a3rsAlternative.a3rs01_0,9);
            Assert.Equal(a3rsExpected.a3rs01_1 ,a3rsAlternative.a3rs01_1,9);
            Assert.Equal(a3rsExpected.a3rs01_2 ,a3rsAlternative.a3rs01_2,9);
            Assert.Equal(a3rsExpected.a3rs02_0 ,a3rsAlternative.a3rs02_0,9);
            Assert.Equal(a3rsExpected.a3rs02_1 ,a3rsAlternative.a3rs02_1,9);
            Assert.Equal(a3rsExpected.a3rs02_2 ,a3rsAlternative.a3rs02_2,9);
            Assert.Equal(a3rsExpected.a3rs10_0 ,a3rsAlternative.a3rs10_0,9);
            Assert.Equal(a3rsExpected.a3rs10_1 ,a3rsAlternative.a3rs10_1,9);
            Assert.Equal(a3rsExpected.a3rs10_2 ,a3rsAlternative.a3rs10_2,9);
            Assert.Equal(a3rsExpected.a3rs11_0 ,a3rsAlternative.a3rs11_0,9);
            Assert.Equal(a3rsExpected.a3rs11_1 ,a3rsAlternative.a3rs11_1,9);
            Assert.Equal(a3rsExpected.a3rs11_2 ,a3rsAlternative.a3rs11_2,9);
            Assert.Equal(a3rsExpected.a3rs12_0 ,a3rsAlternative.a3rs12_0,9);
            Assert.Equal(a3rsExpected.a3rs12_1 ,a3rsAlternative.a3rs12_1,9);
            Assert.Equal(a3rsExpected.a3rs12_2 ,a3rsAlternative.a3rs12_2,9);
            Assert.Equal(a3rsExpected.a3rs20_0 ,a3rsAlternative.a3rs20_0,9);
            Assert.Equal(a3rsExpected.a3rs20_1 ,a3rsAlternative.a3rs20_1,9);
            Assert.Equal(a3rsExpected.a3rs20_2 ,a3rsAlternative.a3rs20_2,9);
            Assert.Equal(a3rsExpected.a3rs21_0 ,a3rsAlternative.a3rs21_0,9);
            Assert.Equal(a3rsExpected.a3rs21_1 ,a3rsAlternative.a3rs21_1,9);
            Assert.Equal(a3rsExpected.a3rs21_2 ,a3rsAlternative.a3rs21_2,9);
            Assert.Equal(a3rsExpected.a3rs22_0 ,a3rsAlternative.a3rs22_0,9);
            Assert.Equal(a3rsExpected.a3rs22_1 ,a3rsAlternative.a3rs22_1,9);
            Assert.Equal(a3rsExpected.a3rs22_2 ,a3rsAlternative.a3rs22_2,9);


            Assert.Equal(da3_drdsExpected[0, 0][0], da3_drds[0, 0][0],9);
            Assert.Equal(da3_drdsExpected[0, 0][1], da3_drds[0, 0][1],9);
            Assert.Equal(da3_drdsExpected[0, 0][2], da3_drds[0, 0][2],9);
            Assert.Equal(da3_drdsExpected[0, 1][0], da3_drds[0, 1][0],9);
            Assert.Equal(da3_drdsExpected[0, 1][1], da3_drds[0, 1][1],9);
            Assert.Equal(da3_drdsExpected[0, 1][2], da3_drds[0, 1][2],9);
            Assert.Equal(da3_drdsExpected[0, 2][0], da3_drds[0, 2][0],9);
            Assert.Equal(da3_drdsExpected[0, 2][1], da3_drds[0, 2][1],9);
            Assert.Equal(da3_drdsExpected[0, 2][2], da3_drds[0, 2][2],9);
                                                                     
            Assert.Equal(da3_drdsExpected[1, 0][0], da3_drds[1, 0][0],9);
            Assert.Equal(da3_drdsExpected[1, 0][1], da3_drds[1, 0][1],9);
            Assert.Equal(da3_drdsExpected[1, 0][2], da3_drds[1, 0][2],9);
            Assert.Equal(da3_drdsExpected[1, 1][0], da3_drds[1, 1][0],9);
            Assert.Equal(da3_drdsExpected[1, 1][1], da3_drds[1, 1][1],9);
            Assert.Equal(da3_drdsExpected[1, 1][2], da3_drds[1, 1][2],9);
            Assert.Equal(da3_drdsExpected[1, 2][0], da3_drds[1, 2][0],9);
            Assert.Equal(da3_drdsExpected[1, 2][1], da3_drds[1, 2][1],9);
            Assert.Equal(da3_drdsExpected[1, 2][2], da3_drds[1, 2][2],9);
                                                                     
            Assert.Equal(da3_drdsExpected[2, 0][0], da3_drds[2, 0][0],9);
            Assert.Equal(da3_drdsExpected[2, 0][1], da3_drds[2, 0][1],9);
            Assert.Equal(da3_drdsExpected[2, 0][2], da3_drds[2, 0][2],9);
            Assert.Equal(da3_drdsExpected[2, 1][0], da3_drds[2, 1][0],9);
            Assert.Equal(da3_drdsExpected[2, 1][1], da3_drds[2, 1][1],9);
            Assert.Equal(da3_drdsExpected[2, 1][2], da3_drds[2, 1][2],9);
            Assert.Equal(da3_drdsExpected[2, 2][0], da3_drds[2, 2][0],9);
            Assert.Equal(da3_drdsExpected[2, 2][1], da3_drds[2, 2][1],9);
            Assert.Equal(da3_drdsExpected[2, 2][2], da3_drds[2, 2][2],9);
            #endregion
        }
    }
}
