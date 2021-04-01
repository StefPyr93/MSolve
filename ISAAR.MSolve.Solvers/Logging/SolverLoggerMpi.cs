using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Distributed;

//TODO: Use enums instead of strings for the solver task and dof category. Or use interfaces & enum classes, to adhere to 
//      open-closed principle.
namespace ISAAR.MSolve.Solvers.Logging
{
    public class SolverLoggerMpi : SolverLoggerBase
    {
        private readonly ProcessDistribution procs;

        public SolverLoggerMpi(ProcessDistribution processDistribution, string solverName) : base(solverName)
        {
            this.procs = processDistribution;
        }

        /// <summary>
        /// Adds the duration of the selected task to the duration of the same task during the current analysis step.
        /// </summary>
        /// <param name="task"></param>
        /// <param name="duration"></param>
        public override void LogCurrentTaskDuration(string task)
        {
            procs.Communicator.Barrier();
            watch.Stop();
            bool exists = taskDurations[CurrentStep].TryGetValue(task, out long durationSofar);
            taskDurations[CurrentStep][task] = durationSofar + watch.ElapsedMilliseconds;
            watch.Reset();
        }
       
        public override void StartMeasuringTime()
        {
            procs.Communicator.Barrier();
            watch.Start();
        }
    }
}
