using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

// Truss nodes:
// 0 -- 1

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationTruss1D
    {
        private static readonly InterpolationTruss1D uniqueInstance = new InterpolationTruss1D();
        private readonly Dictionary<IQuadrature1D, IReadOnlyList<Matrix>> cachedNaturalGradientsAtGPs;
        private readonly Dictionary<IQuadrature1D, IReadOnlyList<double[]>> cachedFunctionsAtGPs;
        private readonly Dictionary<IQuadrature1D, IReadOnlyList<Matrix>> cachedN3AtGPs;

        private InterpolationTruss1D()
        {
            NodalNaturalCoordinates = new NaturalPoint[]
            {
                new NaturalPoint(-1.0),
                new NaturalPoint(+1.0)
            };
            this.cachedNaturalGradientsAtGPs = new Dictionary<IQuadrature1D, IReadOnlyList<Matrix>>();
            this.cachedFunctionsAtGPs = new Dictionary<IQuadrature1D, IReadOnlyList<double[]>>();
            this.cachedN3AtGPs = new Dictionary<IQuadrature1D, IReadOnlyList<Matrix>>();
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public IReadOnlyList<NaturalPoint> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationTruss1D"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationTruss1D UniqueInstance => uniqueInstance;

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        //public override IInverseInterpolation1D CreateInverseMappingFor(IReadOnlyList<Node> nodes)
        //    => new InverseInterpolationTruss1D(nodes);

        public double[] EvaluateAt(double xi)
        {
            var values = new double[2];
            values[0] = 0.50 * (1 - xi);
            values[1] = 0.50 * (1 + xi);
            return values;
        }

        public Matrix EvaluateGradientsAt()
        {
            var derivatives = Matrix.CreateZero(1, 2);
            derivatives[0, 0] = -0.50; // N1,ksi
            derivatives[0, 1] = +0.50; // N2,ksi
            return derivatives;
        }

        public IReadOnlyList<EvalInterpolation1D> EvaluateAllAtGaussPoints(IReadOnlyList<Node> nodes, IQuadrature1D quadrature)
        {
            // The shape functions and natural derivatives at each Gauss point are probably cached from previous calls
            IReadOnlyList<double[]> shapeFunctionsAtGPs = EvaluateFunctionsAtGaussPoints(quadrature);
            IReadOnlyList<Matrix> naturalShapeDerivativesAtGPs = EvaluateNaturalGradientsAtGaussPoints(quadrature);

            // Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
            int numGPs = quadrature.IntegrationPoints.Count;
            var interpolationsAtGPs = new EvalInterpolation1D[numGPs];
            //for (int gp = 0; gp < numGPs; ++gp)
            //{
            //    interpolationsAtGPs[gp] = new EvalInterpolation2D(shapeFunctionsAtGPs[gp],
            //        naturalShapeDerivativesAtGPs[gp], new IsoparametricJacobian2D(nodes, naturalShapeDerivativesAtGPs[gp]));
            //}
            return interpolationsAtGPs;
        }

        private IReadOnlyList<Matrix> EvaluateNaturalGradientsAtGaussPoints(IQuadrature1D quadrature)
        {
            bool isCached = cachedNaturalGradientsAtGPs.TryGetValue(quadrature,
                out IReadOnlyList<Matrix> naturalGradientsAtGPs);
            if (isCached) return naturalGradientsAtGPs;
            else
            {
                int numGPs = quadrature.IntegrationPoints.Count;
                var naturalGradientsAtGPsArray = new Matrix[numGPs];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    GaussPoint gaussPoint = quadrature.IntegrationPoints[gp];
                    naturalGradientsAtGPsArray[gp] = EvaluateGradientsAt();
                }
                cachedNaturalGradientsAtGPs.Add(quadrature, naturalGradientsAtGPsArray);
                return naturalGradientsAtGPsArray;
            }
        }

        public IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature1D quadrature)
        {
            bool isCached = cachedFunctionsAtGPs.TryGetValue(quadrature,
                out IReadOnlyList<double[]> shapeFunctionsAtGPs);
            if (isCached) return shapeFunctionsAtGPs;
            else
            {
                int numGPs = quadrature.IntegrationPoints.Count;
                var shapeFunctionsAtGPsArray = new double[numGPs][];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    GaussPoint gaussPoint = quadrature.IntegrationPoints[gp];
                    shapeFunctionsAtGPsArray[gp] = EvaluateAt(gaussPoint.Xi);
                }
                cachedFunctionsAtGPs.Add(quadrature, shapeFunctionsAtGPsArray);
                return shapeFunctionsAtGPsArray;
            }
        }

        public IReadOnlyList<Matrix> EvaluateN3ShapeFunctionsReorganized(IQuadrature1D quadrature)
        {
            bool isCached = cachedN3AtGPs.TryGetValue(quadrature,
                out IReadOnlyList<Matrix> N3AtGPs);
            if (isCached) return N3AtGPs;
            else
            {
                IReadOnlyList<double[]> N1 = EvaluateFunctionsAtGaussPoints(quadrature);
                N3AtGPs = GetN3ShapeFunctionsReorganized(quadrature, N1);
                cachedN3AtGPs.Add(quadrature, N3AtGPs);
                return N3AtGPs;
            }
        }

        private IReadOnlyList<Matrix> GetN3ShapeFunctionsReorganized(IQuadrature1D quadrature, IReadOnlyList<double[]> N1)
        {
            //TODO reorganize cohesive shell  to use only N1 (not reorganised)

            int nGaussPoints = quadrature.IntegrationPoints.Count;
            var N3 = new Matrix[nGaussPoints]; // shapeFunctionsgpData
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                double ksi = quadrature.IntegrationPoints[npoint].Xi;
                var N3gp = Matrix.CreateZero(3, 6); //8=nShapeFunctions;
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 2; m++) N3gp[l, l + 3 * m] = N1[npoint][m];
                }
                N3[npoint] = N3gp;
            }
            return N3;
        }

        public class EvalInterpolation1D
        {
        }
    }

    
}
