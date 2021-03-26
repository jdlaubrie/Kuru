
import numpy as np
cimport numpy as np

ctypedef double Real
ctypedef long long int Int8 

cdef extern from "_dMUMPS_.h":
    void dMUMPS_job6(Real *rhs, Real *a, int *irn, int *jcn, const int n, const Int8 nnz) nogil

def dMUMPS_solve(np.ndarray[Real, ndim=1, mode='c'] a,
    np.ndarray[int, ndim=1, mode='c'] irn, np.ndarray[int, ndim=1, mode='c'] jcn,
    np.ndarray[Real, ndim=1, mode='c'] rhs):
    
    """ Define the dimensions (n,nnz) and matrix coordinates irn,jcn 
    n:number of rows, nnz:number of nonzero in matrix A, 
    irn:row position, jcn:column position on matrix A """
  
    cdef int n = rhs.shape[0]
    cdef Int8 nnz = a.shape[0]
    
    # MUMPS starts numeration from 1
    irn += 1
    jcn += 1
    
    dMUMPS_job6(&rhs[0], &a[0], &irn[0], &jcn[0], n, nnz)

    return rhs
