using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface ICohesiveZoneMaterial3D : IFiniteElementMaterial
    {
        double[] Tractions { get; }
        IMatrixView ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);

        void ClearTractions();
        new ICohesiveZoneMaterial3D Clone();
    }
}
