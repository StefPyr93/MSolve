using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ThermalMaterial : IFiniteElementMaterial
    {
        public ThermalMaterial(double density, double specialHeatCoeff, double thermalConductivity)
        {
            this.Density = density;
            this.SpecialHeatCoeff = specialHeatCoeff;
            this.ThermalConductivity = thermalConductivity;
        }

        public double Density { get; }
        public double SpecialHeatCoeff { get; }
        public double ThermalConductivity { get; }

        public int ID => throw new NotImplementedException();

        public bool Modified => false;

        public double[] Coordinates { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public double YoungModulus => throw new NotImplementedException(); //TODO: These 2 need redesign urgently.
        public double PoissonRatio => throw new NotImplementedException();

        public void ClearState() { }

        public void ClearStresses() { }

        object ICloneable.Clone() => Clone();
        public ThermalMaterial Clone() => new ThermalMaterial(Density, SpecialHeatCoeff, ThermalConductivity);

        public void ResetModified() { }
        public void SaveState() { }
    }
}
