#from __future__ import division
import numpy as np
from sys import exit
from numpy import einsum
from .MaterialBase import Material
from Kuru.Tensor import trace, Voigt #, makezero


class ArterialWallMixture(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/2*k1m/k2m*(exp(k2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/2*k1c/k2c*(exp(k2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa/2*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where could be possible
        more than one fiber family.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixture, self).__init__(mtype, ndim, **kwargs)

        self.nvar = self.ndim
        self.is_transversely_isotropic = True
        self.energy_type = "internal_energy"
        self.nature = "nonlinear"
        self.fields = "mechanics"

        if self.ndim==3:
            self.H_VoigtSize = 6
        else:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        self.has_low_level_dispatcher = True
        #self.has_low_level_dispatcher = False

        # FIELD VARIABLES AS GROWTH_&_REMODELING AND/OR DEPOSITION STRETCHES, ETC
        self.has_growth_remodeling = True
        self.has_state_variables = True

    def KineticMeasures(self,F, elem=0):
        # first three direction are use for rotation matrix (Normal, Tangential and Axial)
        anisotropic_orientations = self.anisotropic_orientations[elem,:,:]
        # State variables are filled first with remodeling then with densities and growth
        # remodeling (0-8 is elastin, 9 muscle, 10-13 collagen), densities (14-19) growth (20) 
        state_variables = self.StateVariables
        from Kuru.MaterialLibrary.LLDispatch._ArterialWallMixture_ import KineticMeasures
        return KineticMeasures(self,np.ascontiguousarray(F),
                np.ascontiguousarray(anisotropic_orientations),
                np.ascontiguousarray(state_variables))

    def _ConstituentMeasures_(self,F,elem=0):
        anisotropic_orientations = self.anisotropic_orientations[elem,:,:]
        StateVariables = self.StateVariables
        from Kuru.MaterialLibrary.LLDispatch._ArterialWallMixture_ import LLConstituentMeasures
        return LLConstituentMeasures(self,np.ascontiguousarray(F),
                np.ascontiguousarray(anisotropic_orientations),
                np.ascontiguousarray(StateVariables))

    def Hessian(self,StrainTensors,elem=0,gcounter=0):

        #print(np.finfo(np.float64).eps)
        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Normal = self.anisotropic_orientations[elem][0][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Tangential = self.anisotropic_orientations[elem][1][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Axial = self.anisotropic_orientations[elem][2][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Growth tensor (same for all constituent)
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.StateVariables[gcounter][20]*outerNormal + outerTangential
        # Remodeling tensor for elastin
        F_r = I + (self.StateVariables[gcounter][0:9].reshape(3,3) - I)*self.factor_increment
        F_r = np.dot(Rotation.T,np.dot(F_r,Rotation))
        # Elastin Growth and Remodeling tensor
        F_gr = np.dot(F_r,F_g)
        F_gr_inv = np.linalg.inv(F_gr)

        # ELASTIN
        F_e = np.dot(F,F_gr_inv)
        J_e = np.linalg.det(F_e)
        b_e = np.dot(F_e,F_e.T)

        if self.ndim == 3:
            trb = trace(b_e)
        elif self.ndim == 2:
            trb = trace(b_e) + 1.

        H_Voigt = 2.*self.mu*J_e**(-2./3.)*(1./9.*trb*einsum('ij,kl',I,I) - \
                1./3.*(einsum('ij,kl',I,b_e) + einsum('ij,kl',b_e,I)) + \
                1./6.*trb*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)) )*self.StateVariables[gcounter][14]/J

        if self.is_nearly_incompressible:
            H_Voigt += self.pressure*J_e*(einsum('ij,kl',I,I) - \
                (einsum('ik,jl',I,I) + einsum('il,jk',I,I)))*self.StateVariables[gcounter][14]/J
        else:
            H_Voigt += self.kappa*J_e*((2.*J_e-1.)*einsum('ij,kl',I,I) - \
                (J_e-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))*self.StateVariables[gcounter][14]/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling stretch in fibre direction
            lambda_r = 1. + (self.StateVariables[gcounter][fibre_i+8] - 1.)*self.factor_increment
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            if fibre_i is 1:
                k1 = self.k1m*self.StateVariables[gcounter][15]/J
                k2 = self.k2m
                # Anisotropic Stiffness for this key
                H_Voigt += 4.*k1*(1.+2.*k2*(innerFN_e-1.)**2)*np.exp(k2*(innerFN_e-1.)**2)*einsum('ij,kl',outerFN_e,outerFN_e)
                # Active stress for SMC
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]/J
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
                H_Voigt += -2.*(s_act/(self.rho0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*einsum('ij,kl',outerFN,outerFN)
            elif fibre_i is not 1:
                k1 = self.k1c if (innerFN_e-1.0)>=0.0 else 0.075*self.k1c
                k1 = k1*self.StateVariables[gcounter][14+fibre_i]/J
                k2 = self.k2c
                # Anisotropic Stiffness for this key
                H_Voigt += 4.*k1*(1.+2.*k2*(innerFN_e-1.)**2)*np.exp(k2*(innerFN_e-1.)**2)*einsum('ij,kl',outerFN_e,outerFN_e)

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Normal = self.anisotropic_orientations[elem][0][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Tangential = self.anisotropic_orientations[elem][1][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Axial = self.anisotropic_orientations[elem][2][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Growth tensor (same for all constituent)
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.StateVariables[gcounter][20]*outerNormal + outerTangential
        # Remodeling tensor for elastin
        F_r = I + (self.StateVariables[gcounter][0:9].reshape(3,3) - I)*self.factor_increment
        F_r = np.dot(Rotation.T,np.dot(F_r,Rotation))
        # Elastin Growth and Remodeling tensor
        F_gr = np.dot(F_r,F_g)
        F_gr_inv = np.linalg.inv(F_gr)

        # ELASTIN
        F_e = np.dot(F,F_gr_inv)
        J_e = np.linalg.det(F_e)
        b_e = np.dot(F_e,F_e.T)

        if self.ndim == 3:
            trb = trace(b_e)
        elif self.ndim == 2:
            trb = trace(b_e) + 1.

        stress = self.mu*J_e**(-2./3.)*(b_e-1./3.*trb*I)*self.StateVariables[gcounter][14]/J

        if self.is_nearly_incompressible:
            stress += J_e*self.pressure*I*self.StateVariables[gcounter][14]/J
        else:
            stress += self.kappa*J_e*(J_e-1.)*I*self.StateVariables[gcounter][14]/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling stretch in fibre direction
            lambda_r = 1. + (self.StateVariables[gcounter][fibre_i+8] -1.)*self.factor_increment
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            if fibre_i is 1:
                k1 = self.k1m*self.StateVariables[gcounter][15]/J
                k2 = self.k2m
                # Anisotropic Stiffness for this key
                stress += 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e
                # Active stress for SMC
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]/J
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
                stress += (s_act/(self.rho0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN
            elif fibre_i is not 1:
                k1 = self.k1c if (innerFN_e-1.0)>=0.0 else 0.075*self.k1c
                k1 = k1*self.StateVariables[gcounter][14+fibre_i]/J
                k2 = self.k2c
                # Anisotropic Stiffness for this key
                stress += 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e

        return stress

    def ConstituentMeasures(self,StrainTensors,elem=0,gcounter=0):

        density0 = self.rho0
        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]

        #SMC AND COLLAGEN FIBRES
        fibre_stress = np.zeros((5),dtype=np.float64)
        softness = np.zeros((5),dtype=np.float64)
        for fibre_i in [1,2,3,4,5]:
            fraction = self.StateVariables[gcounter][14+fibre_i]/(J*density0)
            if fraction == 0.:
                continue
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling along the fibre
            lambda_r = self.StateVariables[gcounter][8+fibre_i]
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            # Elastic deformation
            innerFN_e = innerFN/lambda_r**2
            if fibre_i is 1:
                k1 = self.k1m*self.StateVariables[gcounter][15]
                k2 = self.k2m
                # Anisotropic Stress for this fibre
                fibre_stress[0] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)
                # Active stress for SMC
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
                fibre_stress[0] += (s_act/density0)*(1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)/(J*fraction)
            elif fibre_i is not 1:
                k1 = self.k1c if (innerFN_e-1.0)>=0.0 else 0.075*self.k1c
                k1 = k1*self.StateVariables[gcounter][14+fibre_i]
                k2 = self.k2c
                # Anisotropic Stress for this fibre
                fibre_stress[fibre_i-1] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)

            # Fibre softness for remodeling
            stiffness = (8.*k2*innerFN_e*(innerFN_e-1.)**2 + 8.*innerFN_e-4.)*k1*\
                np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)/(J*fraction)
            softness[fibre_i-1] = np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness

