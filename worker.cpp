#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <unistd.h>     /* gethostname */
#include "SchurSolve.hpp"

using namespace std;

int main(int argc, char* argv[])
{

    mpi_check(MPI_Init(&argc, &argv));

    //Get cmd arguments
    if(argc != 3)
    {
        cout << "[worker] Usage: " << argv[0] << " N NS"  << endl;
        return 1;
    }
    int N  = int( strtol(argv[1], NULL, 0) );
    int NS = int( strtol(argv[2], NULL, 0) );

    // Show hostname
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Process with PID %d on %s ready to run\n", getpid(), hostname);
    fflush(stdout);

    int rank, mpi_size, err;  
    err = MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    mpi_check(err);
    err = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    mpi_check(err);

    int mpi_parent_size; 
    MPI_Comm parent_comm; 

    //Get inter-communicator to parent process
    err = MPI_Comm_get_parent(&parent_comm); 
    mpi_check(err);
    if (parent_comm == MPI_COMM_NULL)
    {
        cerr << "No parent!" << endl;
        exit(1);
    } 
    err = MPI_Comm_remote_size(parent_comm, &mpi_parent_size); 
    mpi_check(err);
    if (mpi_parent_size != 1)
    {
        cerr << "Something's wrong with the parent" << endl;
        exit(1);
    } 

    /*   
     * The manager is represented as the process with rank 0 in (the remote 
     * group of) parent_comm.  If the workers need to communicate among 
     * themselves, they can use MPI_COMM_WORLD.  
     * 
     * We however merge the parent and child communicators 
     * so that there is only one communicator containing
     * parend and all the worker processes. We will put    
     * parent in front of the workers, so he will have 
     * rank 0.
     */    

    MPI_Comm my_world_comm;
    err = MPI_Intercomm_merge(parent_comm, 1, &my_world_comm);
    mpi_check(err);

    int rank_new, mpi_size_new;  
    err = MPI_Comm_size(my_world_comm, &mpi_size_new);
    mpi_check(err);
    err = MPI_Comm_rank(my_world_comm, &rank_new);
    mpi_check(err);

    cout << "[worker" << rank << "]new rank is: " << rank << "->"<< rank_new
        << " and size is " << mpi_size << "->" << mpi_size_new << endl;


    // child process waits for master to initiate the solution phase
    int pardiso_mtype = -2; // symmetric H_i
    int schur_factorization = 1; //augmented factorization
    int nrhs = 1;
    SchurSolve schurSolver = SchurSolve(pardiso_mtype, schur_factorization, my_world_comm);
    schurSolver.initSystem_OptimalControl(NULL, N, NS);
    schurSolver.solveSystem(NULL, NULL, nrhs);

    MPI_Comm_disconnect(&parent_comm);
    MPI_Comm_free(&my_world_comm);
    MPI_Finalize(); 

    return 0;
}
