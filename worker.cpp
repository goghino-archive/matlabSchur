#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <unistd.h>     /* gethostname */
#include "SchurSolve.hpp"
#include "engine.h"

using namespace std;

int main(int argv, char* argc[])
{

    mpi_check(MPI_Init(&argv, &argc));

    /* workaround to attach GDB */
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Process with PID %d on %s ready to run\n", getpid(), hostname);
    fflush(stdout);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_size < 2) {
        printf("Run with minimum of two processes: mpirun -np 2 [...].\n");
        return -3;
    }


    int N = -1;
    int NS = -1;
    if(argv == 3)
    {
        N = atoi(argc[1]);
        NS = atoi(argc[2]);
    } 

    // if rank == MASTER
    if(mpi_rank == 0)
    {

        if (argv != 3) {
            cout << "Usage: $mpirun -n NP ./solve_problem N NS \n";
            exit(1);
        }

        if (N <= 0 || NS <= 0) {
            printf("Given problem size is invalid.\n");
            exit(2);
        }
   
        //call matlab
        // char comm[500] = "matlab -nosplash -nodisplay -nojvm -nodesktop -r \"mexsolve\"";
        // cout << "Calling the command: " << comm  << " at the MASTER node"<< endl;
        // system(comm);

        Engine* ep;
        if (!(ep = engOpen("")))
        {
            fprintf(stderr, "\n *** Can't start MATLAB engine! *** \n");
            return EXIT_FAILURE;
        }

        engEvalString(ep, "mexsolve");
    }
    else
    {
          // child process waits for master to initiate the solution phase
          int pardiso_mtype = -2; // symmetric H_i
          int schur_factorization = 1; //augmented factorization
          SchurSolve schurSolver = SchurSolve(pardiso_mtype, schur_factorization);
          schurSolver.initSystem_OptimalControl(NULL, N, NS, NULL);
          
          while(1) {
              //do test on termination, set by master process
              int terminate = 0;
              MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD);
              if (terminate)
                  break;
              
              //get flag new_matrix to child processes
              bool new_matrix;
              MPI_Bcast(&new_matrix, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

              //update system if necessary
              if(new_matrix)
              {
                  schurSolver.updateSystem(NULL);
              }

              int nrhs = 1;
              schurSolver.solveSystem(NULL, NULL, nrhs);
              }
    }


    mpi_check(MPI_Finalize());

    return 0;
    }
