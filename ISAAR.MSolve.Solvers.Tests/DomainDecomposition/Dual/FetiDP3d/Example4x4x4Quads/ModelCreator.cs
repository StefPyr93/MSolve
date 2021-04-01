using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads
{
    public static class ModelCreator
    {
        public static Model CreateModel()
        {
            Model model = new Model();
            for (int s = 0; s < 8; ++s)
            {
                model.SubdomainsDictionary[s] = new Subdomain(s);
            }

            (int[] ElementIds, int[] subdomainIds, int[] NodeIds, int[] constraintIds, int[,] ElementNodes, double[,] NodeCoordinates, int[,] SubdElements) =
                GetModelCreationData();

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
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(ElementNodes[i1, j], model.NodesDictionary[ElementNodes[i1, j]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
            }

            for (int i1 = 0; i1 < subdomainIds.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < SubdElements.GetLength(1); i2++)
                {
                    int subdomainID = subdomainIds[i1];
                    Element element = model.ElementsDictionary[SubdElements[i1, i2]];
                    model.SubdomainsDictionary[subdomainID].Elements.Add(element.ID, element);
                }
            }

            for (int i1 = 0; i1 < constraintIds.GetLength(0); i1++)
            {
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

            }

            // Load
            model.Loads.Add(new Load() { Node = model.NodesDictionary[63], DOF = StructuralDof.TranslationZ, Amount = 1.0 });

            return model;
        }

        public static UsedDefinedCornerNodes DefineCornerNodeSelectionSerial(IModel model)
            => new UsedDefinedCornerNodes(DefineCornerNodesSubdomainsAll(model));

        public static Dictionary<ISubdomain, HashSet<INode>> DefineCornerNodesSubdomainsAll(IModel model)
        {
            var cornerNodes = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                cornerNodes[subdomain] = new HashSet<INode>(new INode[] { model.GetNode(63) });
            }
            return cornerNodes;
        }

        public static HashSet<INode> DefineCornerNodesSubdomain(ISubdomain subdomain)
        {
            Debug.Assert(subdomain.ID >= 0 && subdomain.ID < 8);
            return new HashSet<INode>(new INode[] { subdomain.GetNode(63), });
        }

        public static UserDefinedMidsideNodes DefineMidsideNodeSelectionSerial(IModel model)
            => new UserDefinedMidsideNodes(DefineMidsideNodesAll(model), new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ } );

        public static Dictionary<ISubdomain, HashSet<INode>> DefineMidsideNodesAll(IModel model)
        {
            int[] midsideNodeIds = { 62, 64, 58, 70, 38, 100 };

            var midsideNodes = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                midsideNodes[subdomain] = new HashSet<INode>();
                foreach (int midsideNodeId in midsideNodeIds)
                {
                    try
                    {
                        midsideNodes[subdomain].Add(subdomain.GetNode(midsideNodeId));
                    }
                    catch (Exception)
                    {}
                }
            }
            return midsideNodes;
        }

        private static (int[], int[], int[], int[], int[,], double[,], int[,]) GetModelCreationData()
        {
            int[] ElementIds = new int[64] {1,
                17,
                33,
                49,
                5,
                21,
                37,
                53,
                9,
                25,
                41,
                57,
                13,
                29,
                45,
                61,
                2,
                18,
                34,
                50,
                6,
                22,
                38,
                54,
                10,
                26,
                42,
                58,
                14,
                30,
                46,
                62,
                3,
                19,
                35,
                51,
                7,
                23,
                39,
                55,
                11,
                27,
                43,
                59,
                15,
                31,
                47,
                63,
                4,
                20,
                36,
                52,
                8,
                24,
                40,
                56,
                12,
                28,
                44,
                60,
                16,
                32,
                48,
                64};

            int[] subdomainIds = new int[8]{0,
                1,
                2,
                3,
                4,
                5,
                6,
                7};

            //int[] subdomainIds = new int[8]{0,
            //    4,
            //    2,
            //    6,
            //    1,
            //    5,
            //    3,
            //    7};

            int[] NodeIds = new int[125]{1,
                26,
                51,
                76,
                77,
                6,
                31,
                56,
                86,
                87,
                11,
                36,
                61,
                96,
                97,
                16,
                41,
                66,
                106,
                108,
                17,
                42,
                67,
                107,
                109,
                2,
                27,
                52,
                78,
                79,
                7,
                32,
                57,
                88,
                89,
                12,
                37,
                62,
                98,
                99,
                18,
                43,
                68,
                110,
                112,
                19,
                44,
                69,
                111,
                113,
                3,
                28,
                53,
                80,
                81,
                8,
                33,
                58,
                90,
                91,
                13,
                38,
                63,
                100,
                101,
                20,
                45,
                70,
                114,
                116,
                21,
                46,
                71,
                115,
                117,
                4,
                29,
                54,
                82,
                84,
                9,
                34,
                59,
                92,
                94,
                14,
                39,
                64,
                102,
                104,
                22,
                47,
                72,
                118,
                122,
                24,
                49,
                74,
                120,
                124,
                5,
                30,
                55,
                83,
                85,
                10,
                35,
                60,
                93,
                95,
                15,
                40,
                65,
                103,
                105,
                23,
                48,
                73,
                119,
                123,
                25,
                50,
                75,
                121,
                125};

            int[] constraintIds = new int[98]{1,
                26,
                51,
                76,
                77,
                6,
                31,
                56,
                86,
                87,
                11,
                36,
                61,
                96,
                97,
                16,
                41,
                66,
                106,
                108,
                17,
                42,
                67,
                107,
                109,
                2,
                27,
                52,
                78,
                79,
                7,
                89,
                12,
                99,
                18,
                112,
                19,
                44,
                69,
                111,
                113,
                3,
                28,
                53,
                80,
                81,
                8,
                91,
                13,
                101,
                20,
                116,
                21,
                46,
                71,
                115,
                117,
                4,
                29,
                54,
                82,
                84,
                9,
                94,
                14,
                104,
                22,
                122,
                24,
                49,
                74,
                120,
                124,
                5,
                30,
                55,
                83,
                85,
                10,
                35,
                60,
                93,
                95,
                15,
                40,
                65,
                103,
                105,
                23,
                48,
                73,
                119,
                123,
                25,
                50,
                75,
                121,
                125};

            int[,] ElementNodes = new int[64, 8]{{32,31,26,27,7,6,1,2},
                {57,56,51,52,32,31,26,27},
                {88,86,76,78,57,56,51,52},
                {89,87,77,79,88,86,76,78},
                {37,36,31,32,12,11,6,7},
                {62,61,56,57,37,36,31,32},
                {98,96,86,88,62,61,56,57},
                {99,97,87,89,98,96,86,88},
                {43,41,36,37,18,16,11,12},
                {68,66,61,62,43,41,36,37},
                {110,106,96,98,68,66,61,62},
                {112,108,97,99,110,106,96,98},
                {44,42,41,43,19,17,16,18},
                {69,67,66,68,44,42,41,43},
                {111,107,106,110,69,67,66,68},
                {113,109,108,112,111,107,106,110},
                {33,32,27,28,8,7,2,3},
                {58,57,52,53,33,32,27,28},
                {90,88,78,80,58,57,52,53},
                {91,89,79,81,90,88,78,80},
                {38,37,32,33,13,12,7,8},
                {63,62,57,58,38,37,32,33},
                {100,98,88,90,63,62,57,58},
                {101,99,89,91,100,98,88,90},
                {45,43,37,38,20,18,12,13},
                {70,68,62,63,45,43,37,38},
                {114,110,98,100,70,68,62,63},
                {116,112,99,101,114,110,98,100},
                {46,44,43,45,21,19,18,20},
                {71,69,68,70,46,44,43,45},
                {115,111,110,114,71,69,68,70},
                {117,113,112,116,115,111,110,114},
                {34,33,28,29,9,8,3,4},
                {59,58,53,54,34,33,28,29},
                {92,90,80,82,59,58,53,54},
                {94,91,81,84,92,90,80,82},
                {39,38,33,34,14,13,8,9},
                {64,63,58,59,39,38,33,34},
                {102,100,90,92,64,63,58,59},
                {104,101,91,94,102,100,90,92},
                {47,45,38,39,22,20,13,14},
                {72,70,63,64,47,45,38,39},
                {118,114,100,102,72,70,63,64},
                {122,116,101,104,118,114,100,102},
                {49,46,45,47,24,21,20,22},
                {74,71,70,72,49,46,45,47},
                {120,115,114,118,74,71,70,72},
                {124,117,116,122,120,115,114,118},
                {35,34,29,30,10,9,4,5},
                {60,59,54,55,35,34,29,30},
                {93,92,82,83,60,59,54,55},
                {95,94,84,85,93,92,82,83},
                {40,39,34,35,15,14,9,10},
                {65,64,59,60,40,39,34,35},
                {103,102,92,93,65,64,59,60},
                {105,104,94,95,103,102,92,93},
                {48,47,39,40,23,22,14,15},
                {73,72,64,65,48,47,39,40},
                {119,118,102,103,73,72,64,65},
                {123,122,104,105,119,118,102,103},
                {50,49,47,48,25,24,22,23},
                {75,74,72,73,50,49,47,48},
                {121,120,118,119,75,74,72,73},
                {125,124,122,123,121,120,118,119}};

            double[,] NodeCoordinates = new double[125, 3]{{-45,-45,-45},
                {-45,-45,-22.5},
                {-45,-45,0},
                {-45,-45,22.5},
                {-45,-45,45},
                {-45,-22.5,-45},
                {-45,-22.5,-22.5},
                {-45,-22.5,0},
                {-45,-22.5,22.5},
                {-45,-22.5,45},
                {-45,0,-45},
                {-45,0,-22.5},
                {-45,0,0},
                {-45,0,22.5},
                {-45,0,45},
                {-45,22.5,-45},
                {-45,22.5,-22.5},
                {-45,22.5,0},
                {-45,22.5,22.5},
                {-45,22.5,45},
                {-45,45,-45},
                {-45,45,-22.5},
                {-45,45,0},
                {-45,45,22.5},
                {-45,45,45},
                {-22.5,-45,-45},
                {-22.5,-45,-22.5},
                {-22.5,-45,0},
                {-22.5,-45,22.5},
                {-22.5,-45,45},
                {-22.5,-22.5,-45},
                {-22.5,-22.5,-22.5},
                {-22.5,-22.5,0},
                {-22.5,-22.5,22.5},
                {-22.5,-22.5,45},
                {-22.5,0,-45},
                {-22.5,0,-22.5},
                {-22.5,0,0},
                {-22.5,0,22.5},
                {-22.5,0,45},
                {-22.5,22.5,-45},
                {-22.5,22.5,-22.5},
                {-22.5,22.5,0},
                {-22.5,22.5,22.5},
                {-22.5,22.5,45},
                {-22.5,45,-45},
                {-22.5,45,-22.5},
                {-22.5,45,0},
                {-22.5,45,22.5},
                {-22.5,45,45},
                {0,-45,-45},
                {0,-45,-22.5},
                {0,-45,0},
                {0,-45,22.5},
                {0,-45,45},
                {0,-22.5,-45},
                {0,-22.5,-22.5},
                {0,-22.5,0},
                {0,-22.5,22.5},
                {0,-22.5,45},
                {0,0,-45},
                {0,0,-22.5},
                {0,0,0},
                {0,0,22.5},
                {0,0,45},
                {0,22.5,-45},
                {0,22.5,-22.5},
                {0,22.5,0},
                {0,22.5,22.5},
                {0,22.5,45},
                {0,45,-45},
                {0,45,-22.5},
                {0,45,0},
                {0,45,22.5},
                {0,45,45},
                {22.5,-45,-45},
                {22.5,-45,-22.5},
                {22.5,-45,0},
                {22.5,-45,22.5},
                {22.5,-45,45},
                {22.5,-22.5,-45},
                {22.5,-22.5,-22.5},
                {22.5,-22.5,0},
                {22.5,-22.5,22.5},
                {22.5,-22.5,45},
                {22.5,0,-45},
                {22.5,0,-22.5},
                {22.5,0,0},
                {22.5,0,22.5},
                {22.5,0,45},
                {22.5,22.5,-45},
                {22.5,22.5,-22.5},
                {22.5,22.5,0},
                {22.5,22.5,22.5},
                {22.5,22.5,45},
                {22.5,45,-45},
                {22.5,45,-22.5},
                {22.5,45,0},
                {22.5,45,22.5},
                {22.5,45,45},
                {45,-45,-45},
                {45,-45,-22.5},
                {45,-45,0},
                {45,-45,22.5},
                {45,-45,45},
                {45,-22.5,-45},
                {45,-22.5,-22.5},
                {45,-22.5,0},
                {45,-22.5,22.5},
                {45,-22.5,45},
                {45,0,-45},
                {45,0,-22.5},
                {45,0,0},
                {45,0,22.5},
                {45,0,45},
                {45,22.5,-45},
                {45,22.5,-22.5},
                {45,22.5,0},
                {45,22.5,22.5},
                {45,22.5,45},
                {45,45,-45},
                {45,45,-22.5},
                {45,45,0},
                {45,45,22.5},
                {45,45,45}};

            int[,] SubdElements = new int[8, 8]{{1,17,5,21,2,18,6,22},
                {33,49,37,53,34,50,38,54},
                {9,25,13,29,10,26,14,30},
                {41,57,45,61,42,58,46,62},
                {3,19,7,23,4,20,8,24},
                {35,51,39,55,36,52,40,56},
                {11,27,15,31,12,28,16,32},
                {43,59,47,63,44,60,48,64}};

            return (ElementIds, subdomainIds, NodeIds, constraintIds, ElementNodes, NodeCoordinates, SubdElements);
        }
    }
}
