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
#include "SchurSolve.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mpi_check(MPI_Init(NULL, NULL));

    //check number of arguments
    if (nrhs != 3)
    {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                 "Three inputs required.");
         return;
    }

    if (nlhs != 1) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                 "One output required.");
         return;
    }

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    assert(mpi_rank == MASTER);
    
    //make sure the N, NS argument is scalar
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

    //get size of the worker pool
    char fA[500];
    int err = mxGetString(prhs[0], fA, 500);
    int N = (int)mxGetScalar(prhs[1]);
    int NS = (int)mxGetScalar(prhs[2]);

    if ( err != 0 ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:getString",
                "Error in mxGetString");
        return;
    }

    // further parameters
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

#ifdef VERBOSE
    cout << "Exact solution and RHS generated randomly." << endl;
#endif

#ifdef DEBUG_rhs
    RHS.save("b.mat");
#endif

    // Solve KKT system
    int pardiso_mtype = -2; // symmetric H_i
    int schur_factorization = 1; // augmented factorization
    SchurSolve schurSolver = SchurSolve(pardiso_mtype, schur_factorization);
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
