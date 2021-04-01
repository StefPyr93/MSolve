using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads
{
    public static class ModelCreatorInput
    {
        public static Model CreateModel()
        {
            (int[] ElementIds, int[] subdomainIds, int[] NodeIds, int[] constraintIds, int[,] ElementNodes, double[,] NodeCoordinates,
                Dictionary<int, int[]> SubdElements, double[] elementStiffnessFactors) =
                GetModelCreationData();

            Model model = new Model();
            for (int s = 0; s < subdomainIds.GetLength(0); ++s)
            {
                model.SubdomainsDictionary[subdomainIds[s]] = new Subdomain(subdomainIds[s]);
            }

            double E_disp = 3.5; //Gpa
            double ni_disp = 0.4;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };

            for (int i1 = 0; i1 < NodeIds.GetLength(0); i1++)
            {
                int nodeID = NodeIds[i1];
                double nodeCoordX = NodeCoordinates[i1, 0];
                double nodeCoordY = NodeCoordinates[i1, 1];
                double nodeCoordZ = NodeCoordinates[i1, 2];

                model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
            }

            for (int i1 = 0; i1 < ElementIds.GetLength(0); i1++)
            {
                Element e1 = new Element()
                {
                    ID = ElementIds[i1],
                    ElementType = new Hexa8NonLinear(new ElasticMaterial3D() { YoungModulus = elementStiffnessFactors[i1] * E_disp, PoissonRatio = ni_disp }, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(ElementNodes[i1, j], model.NodesDictionary[ElementNodes[i1, j]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
            }

            //for (int i1 = 0; i1 < subdomainIds.GetLength(0); i1++)
            //{
            //    for (int i2 = 0; i2 < SubdElements.GetLength(1); i2++)
            //    {
            //        int subdomainID = subdomainIds[i1];
            //        Element element = model.ElementsDictionary[SubdElements[i1, i2]];
            //        model.SubdomainsDictionary[subdomainID].Elements.Add(element.ID, element);
            //    }
            //}

            foreach (var subdId in SubdElements.Keys)
            {
                foreach(int elementId in SubdElements[subdId])
                {
                    Element element = model.ElementsDictionary[elementId];
                    model.SubdomainsDictionary[subdId].Elements.Add(element.ID, element);
                }
            }

            for (int i1 = 0; i1 < constraintIds.GetLength(0); i1++)
            {
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

            }

            // Load
            var CornerNodesIdAndsubdomains = GetCornerNodesForModel();
            model.Loads.Add(new Load() { Node = model.NodesDictionary[CornerNodesIdAndsubdomains.ElementAt(0).Key], DOF = StructuralDof.TranslationZ, Amount = 1.0 });

            return model;
        }
                
        public static UsedDefinedCornerNodes DefineCornerNodeSelectionSerial(Model model)
            => new UsedDefinedCornerNodes(DefineCornerNodesSubdomainsAll(model));

        public static Dictionary<ISubdomain, HashSet<INode>> DefineCornerNodesSubdomainsAll(Model model)
        {
            var CornerNodesIdAndsubdomains = GetCornerNodesForModel();
            Dictionary<int, HashSet<INode>> cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(CornerNodesIdAndsubdomains,model);

            var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x => x.Item1, x => x.Value);

            //var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);

            return cornerNodes_;
        }


        public static UserDefinedMidsideNodes DefineMidsideNodeSelectionSerial(IModel model)
            => new UserDefinedMidsideNodes(DefineMidsideNodesAll(model), new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ } );

        public static Dictionary<ISubdomain, HashSet<INode>> DefineMidsideNodesAll(IModel model)
        {
            var subdomainOutputPath = (new CnstValues()).exampleOutputPathGen;

            int[] subdExtraConstrsFirstNode = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\subdomain_matrices_and_data\subdExtraConstrsFirstNode.txt");

            Dictionary<int, int[]> CornerNodesAndSubdoaminIds = SamplesConsole.SupportiveClasses.PrintUtilities.ConvertArrayToDictionary(subdExtraConstrsFirstNode);

            var SubdomainsAndMidsideNods = CornerNodesAndSubdoaminIds.Select(x => new KeyValuePair<ISubdomain, HashSet<INode>>( model.GetSubdomain(x.Key), x.Value.Select(y => model.GetNode(y)).ToHashSet())).ToDictionary(x => x.Key, x => x.Value);

            return SubdomainsAndMidsideNods;
        }

        private static (int[], int[], int[], int[], int[,], double[,], Dictionary<int, int[]>,double[]) GetModelCreationData()
        {
            var subdomainOutputPath = (new CnstValues()).exampleOutputPathGen;
            int[] ElementIds = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementIds.txt");

            int[] subdomainIds = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\subdomainIds.txt");
                        
            int[] NodeIds = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeIds.txt");

            int[] constraintIds = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\constraintIds.txt");

            var matReader = new FullMatrixReader(false);

            var elementNodes = matReader.ReadFile(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementNodes.txt");

            int[,] ElementNodes = new int[elementNodes.NumRows, elementNodes.NumColumns];

            for (int i1 = 0; i1 < elementNodes.NumRows; i1++)
            {
                for (int i2 = 0; i2 < elementNodes.NumColumns; i2++)
                {
                    ElementNodes[i1, i2] = (int)elementNodes[i1, i2];
                }
            }

            double[,] NodeCoordinates = matReader.ReadFile(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\NodeCoordinates.txt").CopyToArray2D();

            int[] subdElementData = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\subdElements.txt");

            Dictionary<int, int[]> SubdElements = SamplesConsole.SupportiveClasses.PrintUtilities.ConvertArrayToDictionary(subdElementData);

            double[] elementStiffnessFactors = SamplesConsole.SupportiveClasses.PrintUtilities.ReadVector(subdomainOutputPath + @"\model_overwrite\MsolveModel\" + @"\ElementStiffnessFactors.txt");

            return (ElementIds, subdomainIds, NodeIds, constraintIds, ElementNodes, NodeCoordinates, SubdElements, elementStiffnessFactors);
        }

        public static   Dictionary<int, int[]>  GetCornerNodesForModel()
        {
            var subdomainOutputPath = (new CnstValues()).exampleOutputPathGen;

            int[] cornerNodeAndSubdoaminData = SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(subdomainOutputPath + @"\subdomain_matrices_and_data\CornerNodesAndSubdIds.txt");

            Dictionary<int, int[]> CornerNodesAndSubdoaminIds = SamplesConsole.SupportiveClasses.PrintUtilities.ConvertArrayToDictionary(cornerNodeAndSubdoaminData);

            return CornerNodesAndSubdoaminIds;
        }

        public static Dictionary<int, HashSet<INode>> DefineCornerNodesPerSubdomainAndOtherwise(Dictionary<int, int[]> CornerNodesIdAndsubdomains, Model model)
        {
            // a copy of this is used in modelCreatorInput class
            Dictionary<int, HashSet<INode>> cornerNodesList = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());
            Dictionary<int, HashSet<INode>> cornerNodes = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());

            foreach (Subdomain subdomain in model.EnumerateSubdomains())
            {
                cornerNodesList.Add(subdomain.ID, new HashSet<INode>());
            }

            foreach (int CornerNodeID in CornerNodesIdAndsubdomains.Keys)
            {
                Node node1 = model.NodesDictionary[CornerNodeID];
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            //foreach (Subdomain subdomain in model.Subdomains)
            //{
            //    cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID]);
            //}

            return cornerNodesList;
        }
    }
}
