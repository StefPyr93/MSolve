using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;

//TODO: Use enums instead of strings for the solver task and dof category. Or use interfaces & enum classes, to adhere to 
//      open-closed principle.
namespace ISAAR.MSolve.Solvers.Logging
{
    public class SolverLoggerSerial : SolverLoggerBase
    {
        public SolverLoggerSerial(string solverName) : base(solverName) { }

        /// <summary>
        /// Adds the duration of the selected task to the duration of the same task during the current analysis step.
        /// </summary>
        /// <param name="task"></param>
        /// <param name="duration"></param>
        public override void LogCurrentTaskDuration(string task)
        {
            watch.Stop();
            bool exists = taskDurations[CurrentStep].TryGetValue(task, out long durationSofar);
            taskDurations[CurrentStep][task] = durationSofar + watch.ElapsedMilliseconds;
            watch.Reset();
        }

        public override void StartMeasuringTime() => watch.Start();
    }
}
