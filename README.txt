1. drosos@pilatus01:~/matlabMEX $ module list
   Currently Loaded Modulefiles:
   1) slurm          2) pilatus        3) matlab/r2016a  4) openmpi/2.0.1 

   kardos@icsmaster01 (master):~/misc/matlabMEX$ module list
   Currently Loaded Modulefiles:
     1) use.own         2) gcc/6.1.0       3) openmpi/2.0.1   4) matlab/R2016a

    OpenMPI was configured as following: --prefix=${INSTALLDIR} --enable-mpi-fortran=all --with-pmi --disable-dlopen
    The important flag when using MEX is "--disable-dlopen" which specifies that  all plugins will be slurped into Open
    MPI's libraries and it will cause that Open MPI will not look for / open any DSOs at run time 
    (https://www.open-mpi.org/faq/?category=building#avoid-dso). To elaborate more on this: Open MPI uses a bunch of
    plugins for its functionality.  When you dlopen libmpi in a private namespace (like Matlab does),
    and then libmpi tries to dlopen its  plugins, the plugins can't find the symbols that they need in the main libmpi 
    library (because they're in a private namespace).

2. Set up MATLAB's library preloading

    Copy matlab .rc script into $HOME
    cp $MATLAB_HOME/bin/.matlab7rc.sh ~

    Set LDPATH_PREFIX='/apps/gcc/gcc-6.1.0/lib64/' to point to location where libstdc++.so.6 is located.
    This prevents matlab to load its version of the library from $MATLAB_HOME/sys/os/glnxa64/ directory.
    (or alternatively use hack with LD_PRELOAD in step 6, which is redundant in this case).

3. Tell MATLAB to use MKL libraries

   During Matlab execution LAPACKE routine is invoked from MKL library but base routine is invoked from Matlab library libmwlapack.so. libmwlapack.so supports 8-byte integers. So call of LAPACKE routine with LP64 interface leads to crash.
   
   export BLAS_VERSION=${MKLROOT}/lib/intel64/libmkl_rt.so
   export LAPACK_VERSION=${MKLROOT}/lib/intel64/libmkl_rt.so

4. $ make 

5. Set parameters for MPI modules

   kardos@icsmaster01:$ cat ~/.openmpi/mca-params.conf
   btl_tcp_if_include = eth0
   btl = tcp,sm,self
   pml = ob1

6. $ salloc -N <nworkers + 1>

7. $ make run

   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kardos/privateapps/openmpi/2.0.1/lib \
   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kardos/lib/pardiso \
   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/intel/mkl/lib/intel64 \
   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kardos/PowerGrid/lib \
   LD_PRELOAD="/usr/lib64/libslurm.so /apps/gcc/gcc-6.1.0/lib64/libstdc++.so.6" \
   BLAS_VERSION=/apps/intel/mkl/lib/intel64/libmkl_rt.so \
   LAPACK_VERSION=/apps/intel/mkl/lib/intel64/libmkl_rt.so \
   MKL_INTERFACE_LAYER=LP64 \
   mpirun -np 1 matlab -nojvm -nodisplay -nosplash -r "interface; exit"

