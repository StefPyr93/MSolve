//using System;
//using System.Collections.Generic;
//using System.Text;
//using ISAAR.MSolve.FEM.Entities;
//using Troschuetz.Random;

//namespace ISAAR.MSolve.MSAnalysis.RveTemplatesPaper
//{
//    /// <summary>
//    /// Authors: Odysseas Kokkinos
//    /// </summary>
//    public class RandomCntGeometryGenerator
//    {
//        private readonly int _numberOfElementsPerCnt;
//        private readonly double _cntLength;
//        private readonly double _elementLength;
//        private readonly int _numberOfCnts;
//        private readonly double _standardDeviation;
//        private readonly double _upperAngleBound;
//        private readonly double _matrixLength;
//        private readonly double _matrixWidth;
//        private readonly double _matrixHeight;

//        public RandomCntGeometryGenerator(int numberOfSimulations,
//            int numberOfElementsPerCnt, double cntLength, int numberOfCnts,
//            double standardDeviation, double upperAngleBound,
//            double matrixLength, double matrixWidth, double matrixHeight)
//        {
//            _numberOfElementsPerCnt = numberOfElementsPerCnt;
//            _cntLength = cntLength;
//            if (_numberOfElementsPerCnt <= 0)
//                throw new ArgumentOutOfRangeException("The number of CNTS must be greater than zero");

//            _elementLength = _cntLength / (double) _numberOfElementsPerCnt;
//            _numberOfCnts = numberOfCnts;
//            _standardDeviation = standardDeviation;
//            _upperAngleBound = Math.Abs(upperAngleBound);
//            _matrixLength = matrixLength;
//            _matrixWidth = matrixWidth;
//            _matrixHeight = matrixHeight;
//        }

//        public (int[] nodeIds, double[,] nodeCoordinates, int[,] elementConnectivity) GenerateCnts()
//        {
//            var random = new TRandom();
//            var numberOfNodesPerCnt = _numberOfElementsPerCnt + 1;
//            var nodeIds = new int[_numberOfCnts * numberOfNodesPerCnt];
//            var elementConnectivity = new int[_numberOfCnts * _numberOfElementsPerCnt, 2];
//            var nodalCoordinates = new double[_numberOfCnts * numberOfNodesPerCnt, 3];

//            var counterNode = 0;
//            for (int indexCnt = 0; indexCnt < _numberOfCnts; indexCnt++)
//            {
//                var iNode0 = indexCnt * numberOfNodesPerCnt;
//                nodalCoordinates[iNode0, 0] = random.ContinuousUniform(0.0, _matrixLength);
//                nodalCoordinates[iNode0, 1] = random.ContinuousUniform(0.0, _matrixHeight);
//                nodalCoordinates[iNode0, 2] = random.ContinuousUniform(0.0, _matrixWidth);
//                nodeIds[counterNode] = counterNode;
//                counterNode++;

//                for (int indexElement = 0; indexElement < _numberOfElementsPerCnt; indexElement++)
//                {
//                    double dx;
//                    double dy;
//                    double dz;
//                    var iNode = indexCnt * numberOfNodesPerCnt + indexElement;
//                    if (indexElement == 0)
//                    {
//                        var thetaXY = random.ContinuousUniform(-Math.PI, Math.PI);
//                        var thetaXZ = random.ContinuousUniform(-Math.PI, Math.PI);
//                        dx = Math.Cos(thetaXY) * Math.Cos(thetaXZ) * _elementLength;
//                        dy = Math.Sin(thetaXY) * Math.Cos(thetaXZ) * _elementLength;
//                        dz = Math.Sin(thetaXZ) * _elementLength;
//                    }
//                    else
//                    {
//                        var thetaXY = random.Normal(0.0, 1) * _standardDeviation;
//                        var thetaXZ = random.Normal(0.0, 1) * _standardDeviation;
//                        if (Math.Abs(thetaXY) > _upperAngleBound)
//                            thetaXY = Math.Sign(_upperAngleBound);

//                        if (Math.Abs(thetaXZ) > _upperAngleBound)
//                            thetaXZ = Math.Sign(_upperAngleBound);

//                        var phiXY = Math.Atan2(
//                            nodalCoordinates[indexElement, 1] - nodalCoordinates[indexElement - 1, 1],
//                            nodalCoordinates[indexElement, 0] - nodalCoordinates[indexElement - 1, 0]);
//                        var phiXZ = Math.Atan2(
//                            nodalCoordinates[indexElement, 2] - nodalCoordinates[indexElement - 1, 2],
//                            Math.Sqrt(
//                                Math.Pow(nodalCoordinates[indexElement, 0] - nodalCoordinates[indexElement - 1, 0], 2) +
//                                Math.Pow(nodalCoordinates[indexElement, 1] - nodalCoordinates[indexElement - 1, 1], 2)));

//                        dx = Math.Cos(thetaXY + phiXY) * Math.Cos(thetaXZ + phiXZ);
//                        dy = Math.Sin(thetaXY + phiXY) * Math.Cos(thetaXZ + phiXZ) * _elementLength;
//                        dz = Math.Sin(thetaXZ + phiXZ) * _elementLength;
//                    }

//                    var xNode = nodalCoordinates[iNode, 0] + dx;
//                    var yNode = nodalCoordinates[iNode, 1] + dy;
//                    var zNode = nodalCoordinates[iNode, 2] + dz;

//                    if (xNode > _matrixLength || xNode < 0.0 || 
//                        yNode > _matrixHeight || yNode < 0.0 || 
//                        zNode > _matrixWidth || zNode < 0.0)
//                    {
//                        indexElement-=1;
//                        continue;
//                    }

//                    var distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
//                    var elementId = indexCnt * _numberOfElementsPerCnt + indexElement;
//                    elementConnectivity[elementId, 0] = counterNode - 1;
//                    elementConnectivity[elementId, 1] = counterNode;
//                    nodalCoordinates[iNode + 1, 0] = xNode;
//                    nodalCoordinates[iNode + 1, 1] = yNode;
//                    nodalCoordinates[iNode + 1, 2] = zNode;
//                    nodeIds[counterNode] = counterNode;
//                    counterNode++;

                    
//                }
//            }


//            return (nodeIds, nodalCoordinates, elementConnectivity);
//        }
//    }
//}