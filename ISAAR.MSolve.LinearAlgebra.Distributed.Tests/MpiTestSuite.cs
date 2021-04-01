using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed.Exceptions;
using MPI;

//TODO: Needs printing the stack trace when a test fails. ASAP.
//TODO: Avoid forcing the user to pass the class and method names. At least find a better way to pass the method name.
//TODO: Perhaps the static method should be moved to Solvers.Tests.Program.Main(string[] args) instead of being here and called
//      by SamplesConsole.Program.Main(string[] args).
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Tests
{
    public class MpiTestSuite
    {
        private readonly List<(Action<int> test, string className, string methodName)> tests =
            new List<(Action<int> test, string className, string methodName)>();

        public void AddFact(Action<int> test)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((test, className, methodName));
        }

        public void AddTheory<TInput>(Action<int, TInput> test, TInput input)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((numProcesses => test(numProcesses, input), className, $"{methodName}(args = {input})"));
        }

        public void AddTheory<TInput0, TInput1>(Action<int, TInput0, TInput1> test, TInput0 input0, TInput1 input1)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((numProcesses => test(numProcesses, input0, input1), className, 
                $"{methodName}(args = {input0}, {input1})"));
        }

        public void AddTheory<TInput0, TInput1, TInput2>(Action<int, TInput0, TInput1, TInput2> test,
            TInput0 input0, TInput1 input1, TInput2 input2)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((numProcesses => test(numProcesses, input0, input1, input2), className,
                $"{methodName}(args = {input0}, {input1}, {input2})"));
        }

        public void AddTheory<TInput0, TInput1, TInput2, TInput3>(Action<int, TInput0, TInput1, TInput2, TInput3> test, 
            TInput0 input0, TInput1 input1, TInput2 input2, TInput3 input3)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((numProcesses => test(numProcesses, input0, input1, input2, input3), className,
                $"{methodName}(args = {input0}, {input1}, {input2}, {input3})"));
        }

        public void AddTheory<TInput0, TInput1, TInput2, TInput3, TInput4>(
            Action<int, TInput0, TInput1, TInput2, TInput3, TInput4> test, 
            TInput0 input0, TInput1 input1, TInput2 input2, TInput3 input3, TInput4 input4)
        {
            string methodName = test.Method.Name;
            string className = test.Method.DeclaringType.Name;
            tests.Add((numProcesses => test(numProcesses, input0, input1, input2, input3, input4), className,
                $"{methodName}(args = {input0}, {input1}, {input2}, {input3}, {input4})"));
        }

        public void RunTests(string[] args)
        {
            using (new MPI.Environment(ref args))
            {
                Intracommunicator comm = Communicator.world;
                int numProcesses = int.Parse(args[0]);
                string header = $"Process {comm.Rank}: ";

                Console.WriteLine(header + "Starting running tests.");
                for (int t = 0; t < tests.Count; ++t)
                {
                    (Action<int> test, string className, string methodName) = tests[t];
                    comm.Barrier();
                    try
                    {
                        test(numProcesses);
                        Console.WriteLine(header + $"Test {t} - {className}.{methodName} PASSED!");
                    }
                    catch (MpiProcessesException ex)
                    {
                        Console.WriteLine(header + $"Test {t} - {className}.{methodName} SKIPPED!"
                            + $" Incorrect number ({numProcesses}) of MPI processes: " + ex.Message);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine(header + $"Test {t} - {className}.{methodName} FAILED! \n" + ex.StackTrace);
                    }
                }
                comm.Barrier();
                Console.WriteLine(header + "All tests were completed.");
            }
        }
    }
}
