/*
   PowerGrid Schur solver.

   Reads block diagonal system resulting from stochastic
   PDE problem and solves it. Right-hand side vector is
   generated randomly (based on random solution vector).

   Diagonal blocks are distributed to child processes
   that compute local contribution to Schur complement.
   These are reduced to master process that forms Schur
   complement and computes its part of solution, that is
   sent to child processes that finalize computation by
   computing solution of their part of the system.

*/

#include <mpi.h>
#include <cassert>
#include <mex.h>   
#include <unistd.h>

#include "SchurSolve.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char name[256];
    gethostname(name, 256);
    mexPrintf("[manager]Runing on node %s\n", name);

    int err = MPI_Init(NULL, NULL);
    mpi_check(err);

    int mpi_rank, mpi_size;
    err = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    mpi_check(err);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    mpi_check(err);

    assert(mpi_rank == MASTER);
    if (mpi_size != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:mpiSize",
                "Please run with a single process!");
        return;
    }

    //check number of arguments
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                "Four inputs required: fA, N, NS, np.");
        return;
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                "One output required.");
        return;
    }

    //make sure the N, NS, np arguments are scalars and fA is string
    if ( !mxIsChar(prhs[0]) ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notString",
                "A must be a string.");
        return;
    }
    if ( !mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                "N must be a scalar.");
        return;
    }
    if ( !mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                "NS must be a scalar.");
        return;
    }
    if ( !mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                "np must be a scalar.");
        return;
    }

    //get cmd parameters
    char fA[2048];
    err = mxGetString(prhs[0], fA, 2048);
    int N = (int)mxGetScalar(prhs[1]);
    int NS = (int)mxGetScalar(prhs[2]);
    int worker_size = (int)mxGetScalar(prhs[3]);

    if ( err != 0 ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:getString",
                "Error in mxGetString");
        return;
    }

    if (worker_size < 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:mpiSize",
                "np must be greater than 0!");
        return;
    }

    /*  
     * Now spawn the workers. Note that there is a run-time determination 
     * of what type of worker to spawn, and presumably this calculation must 
     * be done at run time and cannot be calculated before starting 
     * the program. If everything is known when the application is  
     * first started, it is generally better to start them all at once 
     * in a single MPI_COMM_WORLD.  
     */ 
    
    MPI_Comm everyone_comm;  //intercommunicator to workers
    const char* worker_program = "./worker"; //name of worker binary
    std::string N_s = std::to_string(N);
    char *N_ch = (char *)N_s.c_str();  //use char const* as type
    std::string NS_s = std::to_string(NS);
    char *NS_ch = (char *)NS_s.c_str();  //use char const* as type
    char *worker_args[3] = {N_ch, NS_ch, NULL};
    //char *worker_args[3];
    //worker_args[0] = "10";
    //worker_args[1] = "2";
    //worker_args[2] = NULL;
    MPI_Info host_info = MPI_INFO_NULL;
    // MPI_Info_create(&host_info); 
    // MPI_Info_set(host_info, "host", "icsnode13,icsnode15");

    err = MPI_Comm_spawn(worker_program, worker_args, worker_size,  
            host_info, 0, MPI_COMM_SELF, &everyone_comm,  
            MPI_ERRCODES_IGNORE);
    mpi_check(err);

    int size_all;
    err = MPI_Comm_remote_size(everyone_comm,&size_all);
    mpi_check(err);

    /* 
     * Parallel code here. The communicator "everyone_comm" can be used 
     * to communicate with the spawned processes, which have ranks 0,.. 
     * MPI_worker_size-1 in the remote group of the intercommunicator 
     * "everyone_comm". 
     */

    /* 
     * We however merge the parent and child communicators 
     * so that there is only one communicator containing
     * parend and all the worker processes. We will put    
     * parent in front of the workers, so he will have 
     * rank 0.
     */    

    MPI_Comm my_world_comm;
    err = MPI_Intercomm_merge(everyone_comm, 0, &my_world_comm);
    mpi_check(err);

    int rank_new, mpi_size_new;  
    err = MPI_Comm_size(my_world_comm, &mpi_size_new);
    mpi_check(err);
    err = MPI_Comm_rank(my_world_comm, &rank_new);
    mpi_check(err);

    mexPrintf("[manager]my rank is %d and size is %d->%d\n", rank_new, mpi_size, mpi_size_new);

    // further parameters for Schur solver
    bool use_direct_solvers   = false;
    int number_of_rhs = 3;

    CSRdouble* KKT = NULL;
    Matrix2D<double> RHS;
    Matrix2D<double> X;
    Matrix2D<double> X_exact;

    int* index_a = NULL;

    // load KKT and generate RHS and reference exact solution
    // load IPOPT KKT matrix from file        
    KKT = new CSRdouble();
    bool hasCols = false;
    KKT->loadFromFile(fA, hasCols);
    KKT->fillSymmetricNew(index_a);

    int nrows =  KKT->nrows;

    // allocate RHS and solution vector
    RHS.allocate(nrows, number_of_rhs, 1.0);
    X_exact.allocate(nrows, number_of_rhs, 1.0);
    X.allocate(nrows, number_of_rhs, 1.0);

    // generate exact solution and rhs vector or read from file
    X_exact.rinit_symmetric();
    for(int i = 0; i < number_of_rhs; i++)
    {
        // RHS = KKT * X_exact
        KKT->multiply(&X_exact.data[nrows*i], &RHS.data[nrows*i]); 
    }

    // Solve KKT system
    int pardiso_mtype = -2; // symmetric H_i
    int schur_factorization = 1; // augmented factorization
    SchurSolve schurSolver = SchurSolve(pardiso_mtype, schur_factorization, my_world_comm);
    schurSolver.initSystem_OptimalControl(KKT, N, NS, index_a);

    // Only master contains RHS with actual data at this point
    // it is communicated to children inside the solve
    schurSolver.solveSystem(X.data, RHS.data, number_of_rhs);
    schurSolver.errorReport(number_of_rhs, *KKT, RHS.data, X.data);
    schurSolver.timingReport();

    //associate outputs
    mxArray *c_out_m = plhs[0] = mxCreateDoubleMatrix(1, KKT->nrows, mxREAL);
    //TODO

    //generate new KKT data, the structure is the same
    // double *new_data;
    // int size;
    // CSRdouble A;
    // bool hasCols = false;
    // A.loadFromFile(fA, hasCols);
    // A.matrixType = SYMMETRIC;
    // int N = A.nonzeros;
    // new_data = new double[N];
    // for (int i = 0; i < N; i++)
    // {
    //     new_data[i] = A.pData[i] + 0;
    //     A.pData[i] = A.pData[i] + 0;
    // }
    // size = N;

    // //update the matrix with new data and solve the system again
    // schurSolver.updateSystem(new_data);
    // schurSolver.solveSystem(X.data, RHS.data, number_of_rhs);
    // schurSolver.errorReport(number_of_rhs, A, RHS.data, X.data);
    // schurSolver.timingReport();
    // delete [] new_data;
    // A.clear();
    // 

    // clean up
    KKT->clear();
    delete KKT;

    // allocated only at master
    RHS.clear();
    X_exact.clear();
    X.clear();

    MPI_Finalize();
}
