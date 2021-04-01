using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Distributed.Exceptions
{
    /// <summary>
    /// The exception that is thrown when the number of MPI processes invoked by the user does not match the number of MPI
    /// processes needed by a method or property.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MpiProcessesException : MpiException
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="MpiProcessesException"/> class.
        /// </summary>
        public MpiProcessesException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="MpiProcessesException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public MpiProcessesException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="MpiProcessesException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public MpiProcessesException(string message, Exception inner) : base(message, inner)
        { }
    }
}
