#cython: profile=False
#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False

from Kuru.Tensor import determinant,daxpy

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport exp, sqrt, pow

ctypedef double Real


cdef extern from "_ArterialWallMixture_.h" nogil:
    cdef cppclass _ArterialWallMixture_[Real]:
        _ArterialWallMixture_() except +
        _ArterialWallMixture_(Real mu,Real kappa,Real k1m,Real k2m,Real k1c,Real k2c,Real rho0,
            Real s_act,Real stretch_m,Real stretch_0,Real stretch_a,Real f_incr) except +
        void SetParameters(Real mu,Real kappa,Real k1m,Real k2m,Real k1c,Real k2c,Real rho0,
            Real s_act,Real stretch_m,Real stretch_0,Real stretch_a,Real f_incr) except +
        void KineticMeasures(Real *Snp, Real* Hnp, int ndim, int ngauss, const Real *Fnp, 
            int nfibre, const Real *Nnp, int nstatv, const Real *SVnp, int near_incomp, const Real pressure) except +


def KineticMeasures(material, np.ndarray[Real,ndim=3,mode='c'] F,
        np.ndarray[Real,ndim=2,mode='c'] anisotropic_orientations,
        np.ndarray[Real,ndim=2,mode='c'] state_variables):

    cdef int ndim = F.shape[2]
    cdef int ngauss = F.shape[0]
    cdef int nfibre = anisotropic_orientations.shape[0]
    cdef int nstatv = state_variables.shape[1]
    cdef np.ndarray[Real, ndim=3, mode='c'] stress, hessian

    stress = np.zeros((ngauss,ndim,ndim),dtype=np.float64)
    if ndim==3:
        hessian = np.zeros((ngauss,6,6),dtype=np.float64)
    elif ndim==2:
        hessian = np.zeros((ngauss,3,3),dtype=np.float64)

    cdef Real pressure = material.pressure
    cdef int near_incomp = int(material.is_nearly_incompressible)

    cdef _ArterialWallMixture_[Real] mat_obj = _ArterialWallMixture_()
    mat_obj.SetParameters(material.mu,material.kappa,material.k1m,material.k2m,material.k1c,
        material.k2c,material.rho0,material.maxi_active_stress,material.maxi_active_stretch,
        material.zero_active_stretch,material.active_stretch,material.factor_increment)
    mat_obj.KineticMeasures(&stress[0,0,0], &hessian[0,0,0], ndim, ngauss, &F[0,0,0], nfibre,
            &anisotropic_orientations[0,0], nstatv, &state_variables[0,0], near_incomp, pressure)

    return stress, hessian


def LLConstituentMeasures(material, np.ndarray[Real,ndim=3,mode='c'] F,
        np.ndarray[Real,ndim=2,mode='c'] anisotropic_orientations,
        np.ndarray[Real,ndim=2,mode='c'] state_variables):

    cdef Real rho0 = material.rho0
    cdef Real k1m = material.k1m
    cdef Real k2m = material.k2m
    cdef Real k1c = material.k1c
    cdef Real k2c = material.k2c
    cdef Real s_act = material.maxi_active_stress
    cdef Real stretch_m = material.maxi_active_stretch
    cdef Real stretch_0 = material.zero_active_stretch
    cdef Real stretch_a = material.active_stretch

    cdef int ngauss = F.shape[0]
    cdef int ndim = F.shape[2]
    cdef int nfibre = anisotropic_orientations.shape[0]
    cdef int nstatv = state_variables.shape[1]

    cdef np.ndarray[Real,ndim=2, mode='c'] fibre_stress = np.zeros((ngauss,nfibre-1),dtype=np.float64)
    cdef np.ndarray[Real,ndim=2, mode='c'] softness = np.zeros((ngauss,nfibre-1),dtype=np.float64)

    ConstituentMeasures_(&fibre_stress[0,0], &softness[0,0], &F[0,0,0], &state_variables[0,0],
        &anisotropic_orientations[0,0], rho0, k1m, k2m, k1c, k2c, s_act, stretch_m, stretch_0, stretch_a,
        ngauss, ndim, nfibre, nstatv)

    return fibre_stress, softness

cdef ConstituentMeasures_(Real *stress, Real *softness, Real *Fnp, Real *SVnp, Real *Nnp, Real rho0,
        Real k1m, Real k2m, Real k1c, Real k2c, Real s_act, Real stretch_m, Real stretch_0, Real stretch_a,
        int ngauss, int ndim, int nfibre, int nstatv):

    cdef int i,j,n,g
    cdef Real J,frac,lambda_r,innerFN,innerFN_e,k1

    cdef Real *FN = <Real*>malloc(sizeof(Real)*ndim)

    for g in range(ngauss):
        # determinant
        J = Fnp[g*ndim*ndim] * Fnp[g*ndim*ndim+ndim+1]*Fnp[g*ndim*ndim+2*ndim+2] - \
            Fnp[g*ndim*ndim] * Fnp[g*ndim*ndim+ndim+2]*Fnp[g*ndim*ndim+2*ndim+1] - \
            Fnp[g*ndim*ndim+1]*Fnp[g*ndim*ndim+ndim] * Fnp[g*ndim*ndim+2*ndim+2] + \
            Fnp[g*ndim*ndim+1]*Fnp[g*ndim*ndim+ndim+2]*Fnp[g*ndim*ndim+2*ndim]   + \
            Fnp[g*ndim*ndim+2]*Fnp[g*ndim*ndim+ndim] * Fnp[g*ndim*ndim+2*ndim+1] - \
            Fnp[g*ndim*ndim+2]*Fnp[g*ndim*ndim+ndim+1]*Fnp[g*ndim*ndim+2*ndim]
        #SMC AND COLLAGEN FIBRES
        for n in range(1,nfibre):
            frac = SVnp[g*nstatv+14+n]/(J*rho0)
            if frac == 0.:
                continue
            # Fibre direction
            for i in range(ndim):
                FN[i] = 0.0
                for j in range(ndim):
                    FN[i] += Fnp[g*ndim*ndim+i*ndim+j]*Nnp[n*ndim+j]
            # Remodeling along the fibre
            lambda_r = SVnp[g*nstatv+8+n]
            # TOTAL deformation
            innerFN = FN[0]*FN[0] + FN[1]*FN[1] + FN[2]*FN[2]
            # Elastic deformation
            innerFN_e = innerFN/pow(lambda_r,2)
            if n == 1:
                # Anisotropic Stress for this fibre
                stress[g*(nfibre-1)] = 2.*k1m*SVnp[g*nstatv+15]*(innerFN_e-1.)*\
                    exp(k2m*pow((innerFN_e-1.),2))*innerFN_e/(J*frac)
                # Active stress for SMC
                stress[g*(nfibre-1)] += (s_act*SVnp[g*nstatv+15]/rho0)*\
                    (1.-pow(((stretch_m-stretch_a)/(stretch_m-stretch_0)),2))/(J*frac)
                # Fibre softness for remodeling
                stiffness = (8.*k2m*innerFN_e*pow((innerFN_e-1.),2) + 8.*innerFN_e-4.)*\
                    k1m*SVnp[g*nstatv+15]*sqrt(innerFN_e)*exp(k2m*pow((innerFN_e-1.),2))/(J*frac)
            elif n != 1:
                k1 = k1c if (innerFN_e-1.0)>=0.0 else 0.075*k1c
                # Anisotropic Stress for this fibre
                stress[g*(nfibre-1)+n-1] = 2.*k1*SVnp[g*nstatv+14+n]*(innerFN_e-1.)*\
                    exp(k2c*pow((innerFN_e-1.),2))*innerFN_e/(J*frac)
                # Fibre softness for remodeling
                stiffness = (8.*k2c*innerFN_e*pow((innerFN_e-1.),2) + 8.*innerFN_e-4.)*\
                    k1*SVnp[g*nstatv+14+n]*sqrt(innerFN_e)*exp(k2c*pow((innerFN_e-1.),2))/(J*frac)

            softness[g*(nfibre-1)+n-1] = sqrt(innerFN)/(innerFN_e*stiffness)

    free(FN)

