/*
 *  This solver function is based on MUMPS c_example.c
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

/* Function targets JOB=6 of dMUMPS 
   (analyse, factorization and solve)*/
void dMUMPS_job6(double *rhs, double *a,
    MUMPS_INT *irn, MUMPS_INT *jcn,
    const MUMPS_INT n, const MUMPS_INT8 nnz)
{
  DMUMPS_STRUC_C id;

/* When compiling with -DINTSIZE64, MUMPS_INT is 64-bit but MPI
   ilp64 versions may still require standard int for C interface. */
/* MUMPS_INT myid, ierr; */
  int myid, ierr;

  int argc; char **argv;
  int error = 0;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  id.comm_fortran=USE_COMM_WORLD;
  id.par=1; id.sym=0;
  id.job=JOB_INIT;
  dmumps_c(&id);

  /* Define the problem on the host */
  if (myid == 0) {
    id.n = n; id.nnz =nnz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;
  }
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
  id.ICNTL(1)=-1; id.ICNTL(2)=-1;
  id.ICNTL(3)=-1; id.ICNTL(4)=0;

  /* Call the MUMPS package (analyse, factorization and solve). */
  id.job=6;
  dmumps_c(&id);
  if (id.infog[0]<0) {
    printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
        myid, id.infog[0], id.infog[1]);
    error = 1;
  }

  /* Terminate instance. */
  id.job=JOB_END;
  dmumps_c(&id);
  if (myid == 0) {
    if (!error) {
      ierr = MPI_Finalize();
      return;
    } else {
      printf("An error has occured, please check error code returned by MUMPS.\n");
    }
  }
  ierr = MPI_Finalize();
}
