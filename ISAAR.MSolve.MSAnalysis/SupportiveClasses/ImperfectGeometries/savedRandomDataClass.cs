using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    public class savedRandomDataClass
    {
        public double[] random_data_from_file { get; set; }
        private int retrieve_counter = 0;

        public savedRandomDataClass()
        {

        }
        public savedRandomDataClass(double[] random_data_from_file)
        {
            this.random_data_from_file = random_data_from_file;
        }

        public double rand()
        {
            retrieve_counter += 1;
            return random_data_from_file[retrieve_counter];
        }
    }
}
