using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Distributed;

namespace ISAAR.MSolve.Analyzers.ObjectManagers
{
    public class MaterialManagerMpi: IMaterialManager
    {
        private IContinuumMaterial3DDefGrad chosenMaterial;
        private ProcessDistribution procs;
        private List<IContinuumMaterial3DDefGrad> pseudoMaterials = new List<IContinuumMaterial3DDefGrad>();

        private double[][] pseudoMaterialStrains;
        private Dictionary<IContinuumMaterial3DDefGrad, double[]> pseudoMaterialStresses;
        private Dictionary<IContinuumMaterial3DDefGrad, IMatrixView> pseudoMaterialConsMatrices;

        private Dictionary<IContinuumMaterial3DDefGrad, int> pseudoMaterialsMappingToDatabase;

        

        private int[] zbased_globalFromLocal;  // numbering
        private IContinuumMaterial3DDefGrad[] localMaterialList;
        private double[][] localMaterialStrains;
        private double[][] localMaterialStresses;
        private IMatrixView[] localMaterialConsMatrices;
        private int numFullProcesses;
        private int numMaterialsOfFullProcesses;
        private int numMaterialsOfLastProcess;

        private int sizeOfGLobalMaterialDatabase;

        

        public void Initialize()
        {
            sizeOfGLobalMaterialDatabase = 0;
            if (procs.IsMasterProcess)
            { sizeOfGLobalMaterialDatabase = pseudoMaterials.Count; }
            procs.Communicator.Broadcast<int>(ref sizeOfGLobalMaterialDatabase, procs.MasterProcess);

            if (procs.IsMasterProcess)
            {
                BuildTheSimpleMapping();
                pseudoMaterialStresses = pseudoMaterials.Select(x => new KeyValuePair<IContinuumMaterial3DDefGrad, double[]>(x, null)).ToDictionary(x => x.Key, x => x.Value);
                pseudoMaterialConsMatrices = pseudoMaterials.Select(x => new KeyValuePair<IContinuumMaterial3DDefGrad, IMatrixView>(x, null)).ToDictionary(x => x.Key, x => x.Value);
            }

            BuildLocalDataStructures(sizeOfGLobalMaterialDatabase);

            // pithanws me if master ean exei domes pou tha xreiazontai gia thn epikoinwnia me ta remoteMaterials Pseudomaterial
            //BuildMaterials();

            pseudoMaterialStrains = new double[pseudoMaterials.Count()][];

            
            
        }

        private void BuildLocalDataStructures(int sizeOfGLobalMaterialDatabase)
        {
            numFullProcesses = procs.GetNumSubdomainsPerProcess().Count()-1;
            numMaterialsOfFullProcesses = sizeOfGLobalMaterialDatabase / numFullProcesses;
            numMaterialsOfLastProcess = sizeOfGLobalMaterialDatabase % numFullProcesses;

            if(procs.OwnRank==numFullProcesses) //ennooume oti vriskomaste sthn teleftaia 
            {
                zbased_globalFromLocal = new int[numMaterialsOfLastProcess];
                int firstElementValue = sizeOfGLobalMaterialDatabase - numMaterialsOfLastProcess ; // logw zerobased

                for (int i1 = 0; i1 < numMaterialsOfLastProcess; i1++)
                {
                    zbased_globalFromLocal[i1] = firstElementValue + i1;
                }
                localMaterialStrains = new double[numMaterialsOfLastProcess][];
                localMaterialStresses = new double[numMaterialsOfLastProcess][];
                localMaterialConsMatrices = new IMatrixView[numMaterialsOfLastProcess];

                localMaterialList = new IContinuumMaterial3DDefGrad[numMaterialsOfLastProcess];
                for (int i1 = 0; i1 < localMaterialList.Length; i1++)
                {
                    localMaterialList[i1] = (IContinuumMaterial3DDefGrad)chosenMaterial.Clone();
                    localMaterialConsMatrices[i1] = localMaterialList[i1].ConstitutiveMatrix;
                }
                 
            }
            else
            {
                zbased_globalFromLocal = new int[numMaterialsOfFullProcesses];
                int firstElementValue = procs.OwnRank * numMaterialsOfFullProcesses;

                for (int i1 = 0; i1 < numMaterialsOfFullProcesses; i1++)
                {
                    zbased_globalFromLocal[i1] = firstElementValue + i1;
                }
                localMaterialStrains = new double[numMaterialsOfFullProcesses][];
                localMaterialStresses = new double[numMaterialsOfFullProcesses][];
                localMaterialConsMatrices = new IMatrixView[numMaterialsOfFullProcesses];

                localMaterialList = new IContinuumMaterial3DDefGrad[numMaterialsOfFullProcesses];
                for (int i1 = 0; i1 < localMaterialList.Length; i1++)
                {
                    localMaterialList[i1] = (IContinuumMaterial3DDefGrad)chosenMaterial.Clone();
                    localMaterialConsMatrices[i1] = localMaterialList[i1].ConstitutiveMatrix;
                }
            }

            IMatrixView[][] gatheredCons = procs.Communicator.Gather(localMaterialConsMatrices, procs.MasterProcess);

            if (procs.IsMasterProcess)
            {
                var counter = 0;
                for (int i1 = 0; i1 < numFullProcesses + 1; i1++)
                {
                    for (int i2 = 0; i2 < gatheredCons[i1].Length; i2++)
                    {
                        pseudoMaterialConsMatrices[pseudoMaterials.ElementAt(counter)] = gatheredCons[i1][i2];
                        counter++;
                    }
                }
            }

        }

        public void UpdateMaterials()
        {
            MpiUtilities.BroadcastArray(procs.Communicator, ref pseudoMaterialStrains, procs.MasterProcess);

            for (int i1 = 0; i1 < zbased_globalFromLocal.Length; i1++)
            {
                localMaterialStrains[i1] = pseudoMaterialStrains[zbased_globalFromLocal[i1]];// TODO: telika isws den xreiazetai afto to endiameso dianusma
                localMaterialList[i1].UpdateMaterial(localMaterialStrains[i1]);
                localMaterialStresses[i1] = localMaterialList[i1].Stresses;
                localMaterialConsMatrices[i1] = localMaterialList[i1].ConstitutiveMatrix;
            }

            double[][][] gatheredStresses = procs.Communicator.Gather(localMaterialStresses, procs.MasterProcess); // TODO: mapping kateftheian apo afta
            IMatrixView[][] gatheredCons = procs.Communicator.Gather(localMaterialConsMatrices, procs.MasterProcess);

            if (procs.IsMasterProcess)
            {
                var counter = 0;
                for (int i1 = 0; i1 < numFullProcesses + 1; i1++)
                {
                    for (int i2 = 0; i2 < gatheredStresses[i1].Length; i2++)
                    {
                        pseudoMaterialStresses[pseudoMaterials.ElementAt(counter)] = gatheredStresses[i1][i2];
                        pseudoMaterialConsMatrices[pseudoMaterials.ElementAt(counter)] = gatheredCons[i1][i2];
                        counter++;
                    }
                }
            }
        }

        private void BuildTheSimpleMapping()
        {
            pseudoMaterialsMappingToDatabase = new Dictionary<IContinuumMaterial3DDefGrad, int>(pseudoMaterials.Count);
            int counter = 0;
            foreach (var material in pseudoMaterials)
            {
                pseudoMaterialsMappingToDatabase[material] = counter; counter++;
            }
        }

        public MaterialManagerMpi(IContinuumMaterial3DDefGrad coosenMaterial, ProcessDistribution procs)
        {
            this.chosenMaterial = coosenMaterial;
            this.procs = procs;
        }

        //public MaterialManagerMpi(ProcessDistribution procs)
        //{
        //    this.procs = procs;
        //}

        public void AddMaterial(IContinuumMaterial3DDefGrad addedMaterial)
        {
            pseudoMaterials.Add(addedMaterial);
        }

        public void UpdateMaterialStrainForRemote(IContinuumMaterial3DDefGrad pseudoMaterial, double[] pseudoMaterialStrain)
        {
            pseudoMaterialStrains[pseudoMaterialsMappingToDatabase[pseudoMaterial]] = pseudoMaterialStrain;
        }

        public double[] GetMaterialStress(IContinuumMaterial3DDefGrad pseudoMaterial)
        {
            return pseudoMaterialStresses[pseudoMaterial];
        }

        public IMatrixView GetMaterialConstitutiveMatrix(IContinuumMaterial3DDefGrad pseudoMaterial)
        {
            return pseudoMaterialConsMatrices[pseudoMaterial];
        }
        

        

        

        public void SaveState()
        {
            foreach (IContinuumMaterial3DDefGrad material in localMaterialList)
            {
                material.SaveState();
            }
        }


    }
}
