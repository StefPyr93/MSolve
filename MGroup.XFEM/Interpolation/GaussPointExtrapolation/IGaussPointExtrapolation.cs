﻿using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using MGroup.XFEM.Integration.Quadratures;

//TODO: Vector2 and Tensor2D instead of double[]. Either that or always use double[], in which case vectors and tensors are 
//      covered by the same method.
namespace MGroup.XFEM.Interpolation.GaussPointExtrapolation
{
    /// <summary>
    /// Calculates extrapolations of scalar, vector and tensor fields from the integration points of a quadrature (integration 
    /// rule). This can be done at any point, but utility methods for directly outputting the extrapolated fields at the nodes
    /// of finite elements are also provided.
    /// </summary>
    public interface IGaussPointExtrapolation
    {
        /// <summary>
        /// The integration rule which defines the integration points used for extrapolating values and defining an auxiliary 
        /// coordinate system.
        /// </summary>
        IQuadrature Quadrature { get; }

        /// <summary>
        /// Calculates a scalar quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="scalarsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the scalar will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        double ExtrapolateScalarFromGaussPoints(IReadOnlyList<double> scalarsAtGaussPoints, double[] naturalPoint);

        /// <summary>
        /// Calculates a scalar quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="scalarsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        IReadOnlyList<double> ExtrapolateScalarFromGaussPointsToNodes(IReadOnlyList<double> scalarsAtGaussPoints, 
            IIsoparametricInterpolation interpolation);

        /// <summary>
        /// Calculates a tensor quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the tensor will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGaussPoints, double[] naturalPoint);

        /// <summary>
        /// Calculates a tensor quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the tensor will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        double[] ExtrapolateTensorFromGaussPoints(IReadOnlyList<double[]> tensorsAtGaussPoints, double[] naturalPoint);

        /// <summary>
        /// Calculates a tensor quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        IReadOnlyList<double[]> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<double[]> tensorsAtGaussPoints,
            IIsoparametricInterpolation interpolation);

        /// <summary>
        /// Calculates a tensor quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGaussPoints,
            IIsoparametricInterpolation interpolation);

        /// <summary>
        /// Calculates a vector quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the tensor will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        double[] ExtrapolateVectorFromGaussPoints(IReadOnlyList<double[]> vectorsAtGaussPoints, double[] naturalPoint);

        /// <summary>
        /// Calculates a vector quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        IReadOnlyList<double[]> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<double[]> vectorsAtGaussPoints,
            IIsoparametricInterpolation interpolation);
    }
}
