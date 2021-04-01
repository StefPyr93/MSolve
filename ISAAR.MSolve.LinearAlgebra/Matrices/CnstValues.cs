using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public enum ParameterSet { stiffCase, stiffLargerRve};
    public class CnstValues
    {
        private enum PC
        {
            Gerasimos, Cluster, Serafeim
        }

        public static int exampleNo { get; set; }
        
        public static bool useInput_forRVE = true;
        public static bool isInputInCode_forRVE = false;

        public static ParameterSet parameterSet = ParameterSet.stiffCase;

        public static bool useCornerNodesInShell = false;
         

        bool runCluster = true;
        PC computer = PC.Cluster;

        public static bool runOnlyHexaModel { get; set; } = false;
        //public string exampleOutputPathGen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example1\input_matlab";

        public string exampleOutputPathGen
        {
            get
            {
                if (computer == PC.Cluster) { return @"C:\Users\cluster\Documents\Large_rves\examples\example" + exampleNo + @"\input_matlab"; }
                else if (computer == PC.Gerasimos) { return @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example" + exampleNo + @"\input_matlab"; }
                else if (computer == PC.Serafeim) { return @"C:\Users\Serafeim\Desktop\embedding\example" + exampleNo + @"\input_matlab"; }
                else throw new Exception();

                //if (runCluster) { return @"C:\Users\cluster\Documents\Large_rves\examples\example" + exampleNo + @"\input_matlab"; }
                //else { return @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example" + exampleNo + @"\input_matlab"; }
            }
        }
        public string solverPath { get { return exampleOutputPathGen + @"\model_overwrite\subdomain_data_solver"; } }

        public string debugString = @"C:\Users\acivi\Documents\notes_elegxoi_2\develp_3D";

        public string exampleDiscrInputPathGen { get { return exampleOutputPathGen + @"\Msolve_input"; } }

        public string interfaceSolverStatsPath { get { return exampleOutputPathGen + @"\Msolve_solution"; } }

        public string rand_data_vec_path
        {
            get
            {
                if (computer == PC.Cluster) { return @"C:\Users\cluster\Documents\Large_rves\files_from_old_paths\rand_data.txt"; }
                else if (computer == PC.Gerasimos) { return @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_random_geometry_Msolve\REF2_50_000_renu_new_multiple_algorithms_check_develop_copy_for_progr_random_direct_in_C\rand_data.txt"; }
                else if (computer == PC.Serafeim) { return @"C:\Users\Serafeim\Desktop\embedding\files_from_old_paths\rand_data.txt"; }
                else throw new Exception();

                //if (runCluster) { return @"C:\Users\cluster\Documents\Large_rves\files_from_old_paths\rand_data.txt"; }
                //else { return @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_random_geometry_Msolve\REF2_50_000_renu_new_multiple_algorithms_check_develop_copy_for_progr_random_direct_in_C\rand_data.txt"; }
            }
        }
        public CnstValues()
        {

        }

        public static bool printInterfaceSolutionStats { get; set; } = true;

        //-------------------------------------------------
        #region newtohn raphson and homogenisation pcg stats
        public static bool printNRstiffnessMatrices { get; set; } = false;
        public static int printNRiterPreconditioner { get; set; } = -1;
        // public static int analyzerLoadingStep { get; set; } = -1;
        public static int analyzerNRIter { get; set; } = 0;
        public static string analyzerInfo { get; set; }

        public static bool analyzerInfoIsSolutionForNRiters { get; set; } = false;
        public string incrementalPcgStatsOutputFileExtention = @"\interfaceSolver_FetiDP_3d_iterations_Per_Increment.txt";

        public static int stressIncrNo { get; set; } = 0;

        public static bool WriteNRRelatedPcgStats { get; set; } = false;

        public static bool printHomogenizationRHSsAndStiffnessmat1 { get; set; } = false;
        public static int rhsCounter { get; set; } = 1;

        public static bool isFirstStiffnessMatrixPrint { get; set; } = true;

        #endregion

        #region supress output
        /// <summary>
        /// Print globalSolutionStats will be prevented as well.
        /// </summary>
        public static void PreventMATLABandTotalOutput()
        {
            printSolver_run_overwrite_data_region = false;
            print_subdomain_data = false;
            WRITESTIFFNESSES = false;
            //even interface solution stats
            printInterfaceSolutionStats = false;

            print_hexa_model = false;
            run_overwrite_data_region = false;
            writeFetiSolutionVectors = false;

            printGlobalSolutionStats = false;
        }

        public static bool print_subdomain_data { get; set; } = true;

        public static bool WRITESTIFFNESSES { get; set; } = true;

        public static bool print_hexa_model { get; set; } = true;
        public static bool run_overwrite_data_region { get; set; } = true;
        public static bool writeFetiSolutionVectors { get; set; } = true;
        
        #endregion
        public static bool printSolver_run_overwrite_data_region { get; set; } = true;

        public static bool printGlobalSolutionStats { get; set; } = true;

        public static bool printPcgMatRhsEtc_AndInterfaceProblemStats { get; set; } = false; // 

        public static bool printPreconditoner { get; set; } = false; //

        public static bool printGlobalSuiteSparseMatrix { get; set; } = true;

        public void WriteToFileStringArray(string[] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {

                writer.Write(array[i]);
                //writer.Write(' ');

                writer.WriteLine();
            }
            writer.Flush();

            writer.Dispose();
        }

        public static void preventOutputFileWrite()
        {
            SaveDefaultValuesForRestore();
            printInterfaceSolutionStats = false;
            printSolver_run_overwrite_data_region = false;
            printGlobalSolutionStats = false;
            printPcgMatRhsEtc_AndInterfaceProblemStats = false;
            printPreconditoner = false;
        }

        private static bool[] defaultBooleanValues;
        private static void SaveDefaultValuesForRestore()
        {
            defaultBooleanValues = new bool[] {
                printInterfaceSolutionStats ,
                printSolver_run_overwrite_data_region,
                printGlobalSolutionStats,
                printPcgMatRhsEtc_AndInterfaceProblemStats,
                printPreconditoner};
        }

        public static void RestoreDefaultBoolValues()
        {
            printInterfaceSolutionStats = defaultBooleanValues[0];
            printSolver_run_overwrite_data_region = defaultBooleanValues[1];
            printGlobalSolutionStats = defaultBooleanValues[2];
            printPcgMatRhsEtc_AndInterfaceProblemStats = defaultBooleanValues[3];
            printPreconditoner = defaultBooleanValues[4];
        }

        #region macro model booleans and output
        public static bool writeFe2MacroscaleSolution=true;
        #endregion

        public static bool useV2FiniteElements { get; set; } = false;
    }
}
