using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;

namespace MGroup.Stochastic.Structural
{
    public class ModelBuilder
    {
        public Model GetModel(RandomVariable randomVariable, IStochasticDomainMapper domainMapper, int iteration)
        {
            var material = new ElasticMaterial()
            {
                ID = 1,
                PoissonRatio = 0.3,
                YoungModulus = randomVariable.Realize(iteration, domainMapper, null)
            };

            var m = new Model();
            m.NodesDictionary.Add(0, new Node(id: 0, x: 0, y:  0, z: 0 ));
            m.NodesDictionary.Add(1, new Node(id: 1, x: 1, y:  0, z: 0 ));
            m.ElementsDictionary.Add(1, new Element()
            {
                ID = 1,
                ElementType = new EulerBeam3D(material)
            });
            m.Loads.Add(new Load() { Amount = 10, DOF = StructuralDof.TranslationX, Node = m.NodesDictionary[1] });

            return m;
        }
    }
}
