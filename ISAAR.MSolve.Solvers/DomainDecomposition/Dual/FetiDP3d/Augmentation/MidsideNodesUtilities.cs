using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation
{
    public class MidsideNodesUtilities
    {
        private enum AxisDirection { X, Y, Z}

        public static INode[] FindMidsidesOfBrick3D(ISubdomain subdomain, double meshTol = 1E-8)
        {
            INode nodeXminYminZmin = null;
            INode nodeXmaxYminZmin = null;
            INode nodeXminYmaxZmin = null;
            INode nodeXmaxYmaxZmin = null;
            INode nodeXminYminZmax = null;
            INode nodeXmaxYminZmax = null;
            INode nodeXminYmaxZmax = null;
            INode nodeXmaxYmaxZmax = null;
            double xmin = Double.MaxValue;
            double xmax = Double.MinValue;
            double ymin = Double.MaxValue;
            double ymax = Double.MinValue;
            double zmin = Double.MaxValue;
            double zmax = Double.MinValue;

            foreach (INode node in subdomain.EnumerateNodes())
            {
                if (node.X < xmin) xmin = node.X;
                if (node.X > xmax) xmax = node.X;
                if (node.Y < ymin) ymin = node.Y;
                if (node.Y > ymax) ymax = node.Y;
                if (node.Z < zmin) zmin = node.Z;
                if (node.Z > zmax) zmax = node.Z;

                if (xmin == node.X && ymin == node.Y && zmin == node.Z) nodeXminYminZmin = node;
                if (xmax == node.X && ymin == node.Y && zmin == node.Z) nodeXmaxYminZmin = node;
                if (xmin == node.X && ymax == node.Y && zmin == node.Z) nodeXminYmaxZmin = node;
                if (xmax == node.X && ymax == node.Y && zmin == node.Z) nodeXmaxYmaxZmin = node;
                if (xmin == node.X && ymin == node.Y && zmax == node.Z) nodeXminYminZmax = node;
                if (xmax == node.X && ymin == node.Y && zmax == node.Z) nodeXmaxYminZmax = node;
                if (xmin == node.X && ymax == node.Y && zmax == node.Z) nodeXminYmaxZmax = node;
                if (xmax == node.X && ymax == node.Y && zmax == node.Z) nodeXmaxYmaxZmax = node;
            }

            var edges = new BrickEdge[12];
            edges[ 0] = new BrickEdge(AxisDirection.X, new double[] { ymin, zmin }, nodeXminYminZmin, nodeXmaxYminZmin, subdomain, meshTol);
            edges[ 1] = new BrickEdge(AxisDirection.X, new double[] { ymax, zmin }, nodeXminYmaxZmin, nodeXmaxYmaxZmin, subdomain, meshTol);
            edges[ 2] = new BrickEdge(AxisDirection.X, new double[] { ymin, zmax }, nodeXminYminZmax, nodeXmaxYminZmax, subdomain, meshTol);
            edges[ 3] = new BrickEdge(AxisDirection.X, new double[] { ymax, zmax }, nodeXminYmaxZmax, nodeXmaxYmaxZmax, subdomain, meshTol);
                   
            edges[ 4] = new BrickEdge(AxisDirection.Y, new double[] { zmin, xmin }, nodeXminYminZmin, nodeXminYmaxZmin, subdomain, meshTol);
            edges[ 5] = new BrickEdge(AxisDirection.Y, new double[] { zmax, xmin }, nodeXminYminZmax, nodeXminYmaxZmax, subdomain, meshTol);
            edges[ 6] = new BrickEdge(AxisDirection.Y, new double[] { zmin, xmax }, nodeXmaxYminZmin, nodeXmaxYmaxZmin, subdomain, meshTol);
            edges[ 7] = new BrickEdge(AxisDirection.Y, new double[] { zmax, xmax }, nodeXmaxYminZmax, nodeXmaxYmaxZmax, subdomain, meshTol);

            edges[ 8] = new BrickEdge(AxisDirection.Z, new double[] { xmin, ymin }, nodeXminYminZmin, nodeXminYminZmax, subdomain, meshTol);
            edges[ 9] = new BrickEdge(AxisDirection.Z, new double[] { xmax, ymin }, nodeXmaxYminZmin, nodeXmaxYminZmax, subdomain, meshTol);
            edges[10] = new BrickEdge(AxisDirection.Z, new double[] { xmin, ymax }, nodeXminYmaxZmin, nodeXminYmaxZmax, subdomain, meshTol);
            edges[11] = new BrickEdge(AxisDirection.Z, new double[] { xmax, ymax }, nodeXmaxYmaxZmin, nodeXmaxYmaxZmax, subdomain, meshTol);

            return edges.Select(edge => edge.FindMidsideNode()).ToArray();
        }
        
        private class BrickEdge
        {
            private readonly double meshTol;

            public BrickEdge(AxisDirection axis, double[] constCoordinates, INode corner0, INode corner1, ISubdomain subdomain, 
                double meshTol) //TODO: constCoordinates could be infered by the corner nodes!
            {
                this.Axis = axis;
                this.Corners = new INode[] { corner0, corner1 };
                this.meshTol = meshTol;

                Func<INode, bool> chooseNode = null;
                if (axis == AxisDirection.X)
                {
                    chooseNode = n => 
                        (Math.Abs(n.Y - constCoordinates[0]) <= meshTol) 
                        && (Math.Abs(n.Z - constCoordinates[1]) <= meshTol);
                }
                if (axis == AxisDirection.Y)
                {
                    chooseNode = n =>
                        (Math.Abs(n.Z - constCoordinates[0]) <= meshTol)
                        && (Math.Abs(n.X - constCoordinates[1]) <= meshTol);
                }
                if (axis == AxisDirection.Z)
                {
                    chooseNode = n =>
                        (Math.Abs(n.X - constCoordinates[0]) <= meshTol)
                        && (Math.Abs(n.Y - constCoordinates[1]) <= meshTol);
                }
                AllNodes = subdomain.EnumerateNodes().Where(chooseNode).ToArray();
            }

            public AxisDirection Axis { get; }
            public INode[] Corners { get; }
            public INode[] AllNodes { get; }

            public INode FindMidsideNode()
            {
                Func<INode, double> findCoord = null;
                if (Axis == AxisDirection.X) findCoord = n => n.X;
                else if (Axis == AxisDirection.Y) findCoord = n => n.Y;
                else if (Axis == AxisDirection.Z) findCoord = n => n.Z;
                else throw new ArgumentException();

                double max = Corners.Max(findCoord);
                double min = Corners.Min(findCoord);
                double target = double.NaN;
                if (AllNodes.Length % 2 == 1)
                {
                    target = min + 0.5 * (max - min);
                }
                else
                {
                    double interval = (max - min) / (AllNodes.Length - 1);
                    target = min + (AllNodes.Length / 2) * interval;
                }
                return AllNodes.First(n => Math.Abs(findCoord(n) - target) < meshTol);
            }
        }
    }
}
