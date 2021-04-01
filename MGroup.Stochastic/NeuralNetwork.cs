using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace MGroup.Stochastic
{
    public class NeuralNetwork
    {
        Matrix[] weightMatrix;
        Vector[] biasVector;
        Matrix tempMatrix;
        Vector tempVector;

        Vector inputMinValues;
        Vector inputMaxValues;
        Vector outputMinValues;
        Vector outputMaxValues;
        Vector inputRange;
        Vector outputRange;
        Matrix scaleJacobianMatrix;

        public double[][,] weights { get; private set; }
        public double[][] biases { get; private set; }
        public string[] activationFunctions { get; private set; }
        public bool preprocessed { get; set; }
        public double[] constParameters { get; set; }

        public NeuralNetwork()
        {
            preprocessed = false;
        }

        public void ExtractNeuralNetworkParametersFromMatlab()
        {
            string workingDirectory = @"C:\Users\stefp\OneDrive\Desktop\Solution Soimiri\NNparameters";

            string WeightsFileName = "weights.txt";
            string BiasesFileName = "biases.txt";
            string ActivationFunctionsFileName = "activationFunctions.txt";

            string fileNameOnlyWeightsFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(WeightsFileName));
            string fileNameOnlyBiasesFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(BiasesFileName));
            string fileNameOnlyActivationFunctionsFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(ActivationFunctionsFileName));

            string extension_1 = Path.GetExtension(WeightsFileName);
            string extension_2 = Path.GetExtension(BiasesFileName);
            string extension_3 = Path.GetExtension(ActivationFunctionsFileName);

            string currentWeightsFileName = string.Format("{0}{1}", fileNameOnlyWeightsFileName, extension_1);
            string currentBiasesFileName = string.Format("{0}{1}", fileNameOnlyBiasesFileName, extension_2);
            string currentActivationFunctionsFileName = string.Format("{0}{1}", fileNameOnlyActivationFunctionsFileName, extension_3);

            int WeightsRows = File.ReadLines(currentWeightsFileName).Count();
            int BiasesRows = File.ReadLines(currentBiasesFileName).Count();
            int ActivationFunctionsRows = File.ReadLines(currentActivationFunctionsFileName).Count();

            weights = new double[ActivationFunctionsRows][,];
            biases = new double[ActivationFunctionsRows][];
            activationFunctions = new string[ActivationFunctionsRows];

            var numLayer = 0; var numRow = 0;
            var totalRows = new int[ActivationFunctionsRows];
            var totalCols = new int[ActivationFunctionsRows];

            using (TextReader reader = File.OpenText(currentWeightsFileName))
            {
                for (int i = 0; i < WeightsRows; i++)
                {
                    string text = reader.ReadLine();
                    if (text != "") { if (numRow == 0) { string[] bits = text.Split(','); totalCols[numLayer] = bits.Length - 1; } numRow++; }
                    else { totalRows[numLayer] = numRow; numRow = 0; numLayer++; }
                }
            }

            using (TextReader reader = File.OpenText(currentWeightsFileName))
            {
                for (numLayer = 0; numLayer < ActivationFunctionsRows; numLayer++)
                {
                    weights[numLayer] = new double[totalRows[numLayer], totalCols[numLayer]];
                    for (numRow = 0; numRow < totalRows[numLayer]; numRow++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        for (int numCol = 0; numCol < totalCols[numLayer]; numCol++)
                        {
                            weights[numLayer][numRow, numCol] = double.Parse(bits[numCol]);
                        }
                    }
                    reader.ReadLine();
                }
            }

            using (TextReader reader = File.OpenText(currentBiasesFileName))
            {
                for (numLayer = 0; numLayer < ActivationFunctionsRows; numLayer++)
                {
                    biases[numLayer] = new double[totalRows[numLayer]];
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    for (numRow = 0; numRow < totalRows[numLayer]; numRow++)
                    {
                        biases[numLayer][numRow] = double.Parse(bits[numRow]);
                    }
                }
            }

            using (TextReader reader = File.OpenText(currentActivationFunctionsFileName))
            {
                for (numLayer = 0; numLayer < ActivationFunctionsRows; numLayer++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    activationFunctions[numLayer] = bits[0];
                }
            }

            weightMatrix = new Matrix[ActivationFunctionsRows];
            biasVector = new Vector[ActivationFunctionsRows];
            for (numLayer = 0; numLayer < ActivationFunctionsRows; numLayer++)
            {
                weightMatrix[numLayer] = Matrix.CreateFromArray(weights[numLayer]);
                biasVector[numLayer] = Vector.CreateFromArray(biases[numLayer]);
            }
        }

        public void InitializePreprocessing()
        {
            preprocessed = true;
            string workingDirectory = @"C:\Users\stefp\OneDrive\Desktop\Solution Soimiri\NNparameters";
            string InputOutputRangesFileName = "inputoutputRanges.txt";
            string fileNameOnlyInputOutputRangesFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(InputOutputRangesFileName));
            string extension_4 = Path.GetExtension(InputOutputRangesFileName);
            string currentInputOutputRangesFileName = string.Format("{0}{1}", fileNameOnlyInputOutputRangesFileName, extension_4);
            int InputOutputRangesRows = File.ReadLines(currentInputOutputRangesFileName).Count();
            var inputOutputRanges = new double[InputOutputRangesRows][];

            using (TextReader reader = File.OpenText(currentInputOutputRangesFileName))
            {
                for (int i = 0; i < InputOutputRangesRows; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    inputOutputRanges[i] = new double[bits.Length - 1];
                    for (int numCol = 0; numCol < bits.Length - 1; numCol++)
                    {
                        inputOutputRanges[i][numCol] = double.Parse(bits[numCol]);
                    }
                }
            }
            inputMinValues = Vector.CreateFromArray(inputOutputRanges[0]);
            inputMaxValues = Vector.CreateFromArray(inputOutputRanges[1]);
            outputMinValues = Vector.CreateFromArray(inputOutputRanges[2]);
            outputMaxValues = Vector.CreateFromArray(inputOutputRanges[3]);
            inputRange = Vector.CreateFromVector(inputMaxValues);
            inputRange.SubtractIntoThis(inputMinValues);
            outputRange = Vector.CreateFromVector(outputMaxValues);
            outputRange.SubtractIntoThis(outputMinValues);

            var inputRatio = inputRange.DoToAllEntries(x => x / 2);
            var outputRatio = outputRange.DoToAllEntries(x => x / 2);
            inputRatio.DoToAllEntriesIntoThis(x => 1 / x);
            scaleJacobianMatrix = outputRatio.TensorProduct(inputRatio);
            scaleJacobianMatrix = scaleJacobianMatrix.GetSubmatrix(0, 6, 0, 6);
        }


        public Vector CalculateNeuralNetworkOutput(double[] input)
        {
            double[] inputCopy = input.ToArray();
            if (constParameters != null)
                inputCopy = inputCopy.Concat(constParameters).ToArray();
            var outputVector = Vector.CreateFromArray(inputCopy);
            if (preprocessed == true)
            {
                outputVector.SubtractIntoThis(inputMinValues);
                outputVector.DivideEntrywiseIntoThis(inputRange);
                outputVector.ScaleIntoThis(2);
                outputVector.DoToAllEntriesIntoThis(x => x - 1);
            }
            for (int i = 0; i < weightMatrix.Length; i++)
            {
                tempVector = weightMatrix[i].Multiply(outputVector);
                tempVector.AxpyIntoThis(biasVector[i], 1);
                outputVector = CalculateActivationFunction(activationFunctions[i]);
            }
            if (preprocessed == true)
            {
                outputVector.DoToAllEntriesIntoThis(x => x + 1);
                outputVector.DoToAllEntriesIntoThis(x => x/2);
                outputVector.MultiplyEntrywiseIntoThis(outputRange);
                outputVector.AddIntoThis(outputMinValues);
            }
            return outputVector;
        }

        public Matrix CalculateNeuralNetworkJacobian(double[] input)
        {
            double[] inputCopy = input.ToArray();
            if (constParameters != null)
                inputCopy = inputCopy.Concat(constParameters).ToArray();
            var inputVector = Vector.CreateFromArray(inputCopy);
            if (preprocessed == true)
            {
                inputVector.SubtractIntoThis(inputMinValues);
                inputVector.DivideEntrywiseIntoThis(inputRange);
                inputVector.ScaleIntoThis(2);
                inputVector.DoToAllEntriesIntoThis(x => x - 1);
            }
            var derivativeMatrix = Matrix.CreateZero(weights[0].GetLength(0), weights[0].GetLength(1));
            var jacobianMatrix = Matrix.CreateIdentity(weights[0].GetLength(1));

            for (int i = 0; i < weightMatrix.Length; i++)
            {
                tempVector = weightMatrix[i].Multiply(inputVector);
                tempVector.AxpyIntoThis(biasVector[i], 1);
                tempMatrix = weightMatrix[i];
                inputVector = CalculateActivationFunction(activationFunctions[i]);
                derivativeMatrix = CalculateActivationFunctionDerivative(activationFunctions[i]);
                jacobianMatrix = jacobianMatrix.MultiplyLeft(derivativeMatrix);
            }
            jacobianMatrix = jacobianMatrix.GetSubmatrix(0, 6, 0, 6);
            if (preprocessed == true)
            {
                jacobianMatrix.MultiplyEntrywiseIntoThis(scaleJacobianMatrix);
            }
            return jacobianMatrix;
        }

        private Vector CalculateActivationFunction(string activationFunction)
        {
            Vector CalculatedFunction = tempVector.Copy();
            if (activationFunction == "purelin")
            { }
            else if (activationFunction == "tansig")
            {
                CalculatedFunction.DoToAllEntriesIntoThis(x => 2 / (Math.Exp(-2 * x) + 1) - 1);
            }
            return CalculatedFunction;
        }

        private Matrix CalculateActivationFunctionDerivative(string activationFunction)
        {
            var CalculatedFunction = tempMatrix;
            if (activationFunction == "purelin")
            { }
            else if (activationFunction == "tansig")
            {
                var tempVectorDer = tempVector.Copy();
                tempVectorDer.DoToAllEntriesIntoThis(x => Math.Exp(-2 * x));
                tempVectorDer.DoToAllEntriesIntoThis(x => 4 * x / ((x + 1) * (x + 1)));
                //var DiagMatrix = Matrix.CreateIdentity(tempVectorDer.Length);
                //DiagMatrix.SetDiagonal(tempVectorDer);
                var DiagMatrix = Matrix.CreateDiagonalFromVector(tempVectorDer);
                CalculatedFunction = CalculatedFunction.MultiplyLeft(DiagMatrix);
            }
            return CalculatedFunction;
        }
    }
}
