using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.XFEM.Materials
{
    public class HomogeneousElasticMaterial2D: IXMaterialField2D
    {
        private HomogeneousElasticMaterial2D(int id, double youngModulus, double equivalentYoungModulus, double poissonRatio,
            double equivalentPoissonRatio, double thickness)
        {
            this.ID = id;

            MaterialUtilities.CheckYoungModulus(youngModulus);
            MaterialUtilities.CheckPoissonRatio(poissonRatio);
            MaterialUtilities.CheckThickness(thickness);

            HomogeneousYoungModulus = youngModulus;
            HomogeneousEquivalentYoungModulus = equivalentYoungModulus;
            HomogeneousPoissonRatio = poissonRatio;
            HomogeneousEquivalentPoissonRatio = equivalentPoissonRatio;
            HomogeneousThickness = thickness;
        }

        public int ID { get; }
        public double HomogeneousYoungModulus { get; }
        public double HomogeneousEquivalentYoungModulus { get; }
        public double HomogeneousPoissonRatio { get; }
        public double HomogeneousEquivalentPoissonRatio { get; }
        public double HomogeneousThickness { get; }

        public static HomogeneousElasticMaterial2D CreateMaterialForPlaneStrain(int id, double youngModulus,
            double poissonRatio)
        {
            double equivalentE = youngModulus / (1.0 - poissonRatio * poissonRatio);
            double equivalentV = poissonRatio / (1.0 - poissonRatio);
            double thickness = 1.0;
            return new HomogeneousElasticMaterial2D(id, youngModulus, equivalentE, poissonRatio, equivalentV, thickness);
        }

        public static HomogeneousElasticMaterial2D CreateMaterialForPlaneStress(int id, double youngModulus,
            double poissonRatio, double thickness)
        {
            return new HomogeneousElasticMaterial2D(id, youngModulus, youngModulus, poissonRatio, poissonRatio, thickness);
        }

        public IXMaterialField2D Clone()
        {
            return new HomogeneousElasticMaterial2D(this.ID, this.HomogeneousYoungModulus, this.HomogeneousEquivalentYoungModulus,
                this.HomogeneousPoissonRatio, this.HomogeneousEquivalentPoissonRatio, this.HomogeneousThickness);
        }

        public double GetYoungModulusAt(NaturalPoint point, EvalInterpolation2D interpolation) => HomogeneousYoungModulus;

        public double GetEquivalentYoungModulusAt(NaturalPoint point, EvalInterpolation2D interpolation)
            => HomogeneousEquivalentYoungModulus;

        public double GetPoissonRatioAt(NaturalPoint point, EvalInterpolation2D interpolation)
            => HomogeneousPoissonRatio;

        public double GetEquivalentPoissonRatioAt(NaturalPoint point, EvalInterpolation2D interpolation)
            => HomogeneousEquivalentPoissonRatio;

        public double GetThicknessAt(NaturalPoint point, EvalInterpolation2D interpolation)
            => HomogeneousThickness;
        
        public Matrix CalculateConstitutiveMatrixAt(NaturalPoint point, EvalInterpolation2D interpolation)
        {
            var matrix = Matrix.CreateZero(3, 3);
            double eqE = HomogeneousEquivalentYoungModulus;
            double eqV = HomogeneousEquivalentPoissonRatio;
            double scalar = eqE / (1 - eqV * eqV);
            matrix[0, 0] = scalar;
            matrix[0, 1] = scalar * eqV;
            matrix[1, 0] = scalar * eqV;
            matrix[1, 1] = scalar;
            matrix[2, 2] = 0.5 * eqE / (1 + eqV);
            return matrix;
        }

        public void UpdateDistributions() { } // Do nothing for homogeneous preperties
        public void UpdateStrains() { } // Do nothing for elastic properties.
    }
}
