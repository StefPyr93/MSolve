using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using Troschuetz.Random;

namespace MGroup.Stochastic
{
    class RandomCntGeometryGenerator
    {
        private readonly int numberOfElementsPerCnt;
        private readonly double cntLength;
        private readonly double elementLength;
        private readonly int numberOfCnts;
        private readonly double standardDeviation;
        private readonly double upperAngleBound;
        private readonly double matrixLength;
        private readonly double matrixWidth;
        private readonly double matrixHeight;

        public bool periodicInclusions = false;

        public RandomCntGeometryGenerator(int numberOfSimulations,
        int numberOfElementsPerCnt, double cntLength, int numberOfCnts,
        double matrixLength, double matrixWidth, double matrixHeight)
        {
            this.numberOfElementsPerCnt = numberOfElementsPerCnt;
            this.cntLength = cntLength;
            if (numberOfElementsPerCnt <= 0)
                throw new ArgumentOutOfRangeException("The number of CNTS must be greater than zero");

            this.elementLength = cntLength / (double)numberOfElementsPerCnt;
            this.numberOfCnts = numberOfCnts;
            this.matrixLength = matrixLength;
            this.matrixWidth = matrixWidth;
            this.matrixHeight = matrixHeight;
        }

        public (int[] nodeIds, double[][] nodeCoordinates, int[,] elementConnectivity) GenerateCnts()
        {
            var random = new TRandom();
            var numberOfNodesPerCnt = numberOfElementsPerCnt + 1;
            var nodeIds = new int[numberOfCnts * numberOfNodesPerCnt];
            var elementConnectivity = new int[numberOfCnts * numberOfElementsPerCnt, 2];
            var nodalCoordinates = new double[numberOfCnts * numberOfNodesPerCnt][];

            var counterNode = 0;
            for (int indexCnt = 0; indexCnt < numberOfCnts; indexCnt++)
            {
                var iNode0 = indexCnt * numberOfNodesPerCnt;
                nodalCoordinates[iNode0] = new double[3];
                nodalCoordinates[iNode0][0] = random.ContinuousUniform(0.0, matrixLength);
                nodalCoordinates[iNode0][1] = random.ContinuousUniform(0.0, matrixHeight);
                nodalCoordinates[iNode0][2] = random.ContinuousUniform(0.0, matrixWidth);
                nodeIds[counterNode] = counterNode;
                counterNode++;

                for (int indexElement = 0; indexElement < numberOfElementsPerCnt; indexElement++)
                {
                    var iNode = indexCnt * numberOfNodesPerCnt + indexElement;
                    //if (indexElement == 0)
                    //{
                        var randVector = new double[3] {random.Normal(0.0, 1.0), random.Normal(0.0, 1.0), random.Normal(0.0, 1.0)};
                        var normValue = Math.Sqrt(randVector[1] * randVector[1] + randVector[2] * randVector[2] + randVector[3] * randVector[3]);
                        randVector[1] = randVector[1] / normValue; randVector[2] = randVector[2] / normValue; randVector[3] = randVector[3] / normValue;
                    //}
                    //else
                    //{
                    //}
                    var xNode = new double(); var yNode = new double(); var zNode = new double();
                    if (periodicInclusions == false)
                    { 
                        xNode = nodalCoordinates[iNode][0] + randVector[1] * cntLength / numberOfElementsPerCnt;
                        yNode = nodalCoordinates[iNode][1] + randVector[2] * cntLength / numberOfElementsPerCnt;
                        zNode = nodalCoordinates[iNode][2] + randVector[3] * cntLength / numberOfElementsPerCnt;
                        if (xNode > matrixLength || xNode < 0.0 ||
                            yNode > matrixHeight || yNode < 0.0 ||
                            zNode > matrixWidth || zNode < 0.0)
                        {
                            indexElement -= 1;
                            continue;
                        }
                    }
                    else
                    {
                        (xNode, yNode, zNode) = CreatePeriodicInclusion(nodalCoordinates[iNode], randVector);
                    }
                    //var distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    var elementId = indexCnt * numberOfElementsPerCnt + indexElement;
                    elementConnectivity[elementId, 0] = counterNode - 1;
                    elementConnectivity[elementId, 1] = counterNode;
                    nodalCoordinates[iNode + 1][0] = xNode;
                    nodalCoordinates[iNode + 1][1] = yNode;
                    nodalCoordinates[iNode + 1][2] = zNode;
                    nodeIds[counterNode] = counterNode;
                    counterNode++;
                }
            }
            return (nodeIds, nodalCoordinates, elementConnectivity);
        }

        private (double xNode, double yNode, double zNode) CreatePeriodicInclusion(double[] startNodeCoordinates,double[] randVector)
        {
            var xNode = startNodeCoordinates[1] + randVector[1] * cntLength / numberOfElementsPerCnt;
            var yNode = startNodeCoordinates[2] + randVector[2] * cntLength / numberOfElementsPerCnt;
            var zNode = startNodeCoordinates[3] + randVector[3] * cntLength / numberOfElementsPerCnt;
            var distFromBoundary = new double[3, 2];
            var auxNodeCoordinates = startNodeCoordinates;
            var auxElementLength = cntLength / numberOfElementsPerCnt;
            var figureElements = new List<double[]>();
            var done = false;
            while (done == false)
            {
                distFromBoundary[0, 0] = (matrixLength / 2 - auxNodeCoordinates[1]) * -1 / ((randVector[1] * auxElementLength) * -1);
                distFromBoundary[0, 1] = (-matrixLength / 2 - auxNodeCoordinates[1]) * 1 / ((randVector[1] * auxElementLength) * 1);
                distFromBoundary[1, 0] = (matrixWidth / 2 - auxNodeCoordinates[2]) * -1 / ((randVector[2] * auxElementLength) * -1);
                distFromBoundary[1, 1] = (-matrixWidth / 2 - auxNodeCoordinates[2]) * 1 / ((randVector[2] * auxElementLength) * 1);
                distFromBoundary[2, 0] = (matrixHeight / 2 - auxNodeCoordinates[3]) * -1 / ((randVector[3] * auxElementLength) * -1);
                distFromBoundary[2, 1] = (-matrixHeight / 2 - auxNodeCoordinates[3]) * 1 / ((randVector[3] * auxElementLength) * 1);
                double minDist = 1000;
                for (int i = 0; i < distFromBoundary.GetLength(1); i++)
                {
                    for (int j = 0; j < distFromBoundary.GetLength(2); j++)
                    {
                        if (distFromBoundary[i, j] < minDist)
                        {
                            minDist = distFromBoundary[i, j];
                        }
                    }
                }
                if (distFromBoundary[0, 0] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] - matrixLength + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else if (distFromBoundary[0, 1] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] + matrixLength + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else if (distFromBoundary[1, 0] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] - matrixWidth + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else if (distFromBoundary[1, 1] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] + matrixWidth + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else if (distFromBoundary[2, 0] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] - matrixHeight + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else if (distFromBoundary[2, 1] == minDist)
                {
                    auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
                    auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
                    auxNodeCoordinates[3] = auxNodeCoordinates[3] + matrixHeight + minDist * randVector[3] * auxElementLength;
                    auxElementLength = auxElementLength * (1 - minDist);
                }
                else
                {
                    done = true;
                    xNode = auxNodeCoordinates[1] + randVector[1] * auxElementLength;
                    yNode = auxNodeCoordinates[2] + randVector[2] * auxElementLength;
                    zNode = auxNodeCoordinates[3] + randVector[3] * auxElementLength;
                }               
            }
            return (xNode, yNode, zNode);
        }
    }
}
