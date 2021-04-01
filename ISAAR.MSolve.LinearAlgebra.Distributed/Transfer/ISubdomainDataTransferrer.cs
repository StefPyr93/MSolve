using System;
using System.Collections.Generic;
using System.Text;

//TODO: various strategies, which should probably be provided as distinct methods or different implementations of the same interface: 
//  1) send each data of each subdomain one after the other (even to the same process), pack data of the same cluster and them to 
//  each process, put data of all subdomains (and all clusters) in a single array and scatter it, etc.
//  2) Pack and unpack the data, using a dedicated class/interface or just delegates
//  3) only send/receive data of some subdomains, use a dedicated object to easily query which subdomains
//  4) use asynchronous send/receive or one thread to wait for transfering, while another goes on packing/unpacking data, with 
//     respect to the total available RAM.
//  5) Control serialization/deserialization and combine it with packing/unpacking
//TODO: If the data is an array of primitives, then serialization/deserialization can be avoided. Perhaps a dedicated class is 
//      needed for that, or even better dedicated methods in this interface
//TODO: Perhaps provide overloads that take IEnumerable<ISubdomain> and function pointers for gathering the data before 
//      transmission and using them after.
//TODO: Also have variants of the gather methods with an Action<int,T> to use each received item as soon as it is available and 
//      then discard it. In that case, Dictionary<int, T> is not returned. These methods are to be used if there is not enough
//      memory to gather all data in one process and afterwards use them.
//TODO: Scatter, Gather operations that entail all subdomains should be provided as extension methods. Transferrer classes need 
//      only provide versions where only some subdomains are used.
namespace ISAAR.MSolve.LinearAlgebra.Distributed.Transfer
{
    //TODO: Not sure that the subdomainID is needed
    public delegate TPacked PackSubdomainData<TRaw, TPacked>(int subdomainID, TRaw originalData);
    public delegate TRaw UnpackSubdomainData<TRaw, TPacked>(int subdomainID, TPacked packedData);
    public delegate int GetArrayLengthOfPackedData<TRaw>(int subdomainID, TRaw originalData);
    public delegate void PackSubdomainDataIntoArray<TRaw, TPacked>(int subdomainID, TRaw originalData, 
        TPacked[] buffer, int offset);
    public delegate TRaw UnpackSubdomainDataFromArray<TRaw, TPacked>(int subdomainID, TPacked[] buffer, int start, int end);

    public interface ISubdomainDataTransferrer
    {
        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, T> GatherFromAllSubdomains<T>(Dictionary<int, T> processSubdomainsData);

        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, T[]> GatherFromAllSubdomains<T>(Dictionary<int, T[]> processSubdomainsData);

        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData);

        /// <summary>
        /// Like <see cref="GatherFromAllSubdomainsPacked{TRaw, TPacked}(Dictionary{int, TRaw}, PackSubdomainData{TRaw, TPacked}, 
        /// UnpackSubdomainData{TRaw, TPacked})"/>, but packing is done using arrays of possible primitive data types to avoid
        /// serialization/deserialization.
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, TRaw> GatherFromAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData,
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData,
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData);

        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, T> GatherFromSomeSubdomains<T>(Dictionary<int, T> processSubdomainsData,
            ActiveSubdomains activeSubdomains);

        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, T[]> GatherFromSomeSubdomains<T>(Dictionary<int, T[]> processSubdomainsData,
            ActiveSubdomains activeSubdomains);

        /// <summary>
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData,
            ActiveSubdomains activeSubdomains);

        /// <summary>
        /// Like <see cref="GatherFromSomeSubdomainsPacked{TRaw, TPacked}(Dictionary{int, TRaw}, PackSubdomainData{TRaw, TPacked},
        /// UnpackSubdomainData{TRaw, TPacked}, ActiveSubdomains)"/>, but packing is done using arrays of possible primitive 
        /// data types to avoid serialization/deserialization.
        /// In master process this method returns a Dictionary with the data for each subdomain. In other processes, it returns 
        /// null.
        /// </summary>
        Dictionary<int, TRaw> GatherFromSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> processSubdomainsData,
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData,
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains);

        Dictionary<int, T> ScatterToAllSubdomains<T>(Dictionary<int, T> allSubdomainsData_master);

        Dictionary<int, T[]> ScatterToAllSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master);

        Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData);

        /// <summary>
        /// Like <see cref="ScatterToAllSubdomainsPacked{TRaw, TPacked}(Dictionary{int, TRaw}, PackSubdomainData{TRaw, TPacked}, 
        /// UnpackSubdomainData{TRaw, TPacked})"/>, but packing is done using arrays of possible primitive data types to avoid
        /// serialization/deserialization.
        /// </summary>
        Dictionary<int, TRaw> ScatterToAllSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData,
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData);

        Dictionary<int, T> ScatterToSomeSubdomains<T>(Dictionary<int, T> allSubdomainsData_master,
            ActiveSubdomains activeSubdomains);

        Dictionary<int, T[]> ScatterToSomeSubdomains<T>(Dictionary<int, T[]> allSubdomainsData_master, 
            ActiveSubdomains activeSubdomains);

        Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master, 
            PackSubdomainData<TRaw, TPacked> packData, UnpackSubdomainData<TRaw, TPacked> unpackData, 
            ActiveSubdomains activeSubdomains);

        /// <summary>
        /// Like <see cref="ScatterToSomeSubdomainsPacked{TRaw, TPacked}(Dictionary{int, TRaw}, PackSubdomainData{TRaw, TPacked},
        /// UnpackSubdomainData{TRaw, TPacked}, ActiveSubdomains)"/>, but packing is done using arrays of possible primitive data
        /// types to avoid serialization/deserialization.
        Dictionary<int, TRaw> ScatterToSomeSubdomainsPacked<TRaw, TPacked>(Dictionary<int, TRaw> allSubdomainsData_master,
            GetArrayLengthOfPackedData<TRaw> getPackedDataLength, PackSubdomainDataIntoArray<TRaw, TPacked> packData, 
            UnpackSubdomainDataFromArray<TRaw, TPacked> unpackData, ActiveSubdomains activeSubdomains);
    }
}
