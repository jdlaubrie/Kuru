import numpy as np
from numpy import einsum, asarray, eye
from .MaterialBase import Material
from Kuru.Tensor import trace, Voigt


#####################################################################################################
                                    # NEARLY INCOMPRESSIBLE NEOHOOKEAN
#####################################################################################################


class NeoHookean_2(Material):
    """Material model for nearly incompressible neo-Hookean with the following internal energy:

        W(C) = mu/2*J**(-2/3)*(C:I-3)     # for isochoric part
        U(J) = k/2*(J-1)**2               # for volumetric part

        """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(NeoHookean_2, self).__init__(mtype, ndim, **kwargs)

        self.is_transversely_isotropic = False
        self.energy_type = "internal_energy"
        self.nature = "nonlinear"
        self.fields = "mechanics"

        if self.ndim==3:
            self.H_VoigtSize = 6
        elif self.ndim==2:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        self.has_low_level_dispatcher = True

    def KineticMeasures(self,F, elem=0):
        from Kuru.MaterialLibrary.LLDispatch._NeoHookean_2_ import KineticMeasures
        return KineticMeasures(self,np.ascontiguousarray(F))


    def Hessian(self, StrainTensors, elem=0, gcounter=0):
        """Hessian split into isochoroic and volumetric parts"""

        I = StrainTensors['I']
        b = StrainTensors['b'][gcounter]
        J = StrainTensors['J'][gcounter]

        mu = self.mu

        # ISOCHORIC
        H_Voigt = 2*mu*J**(-5./3.)*(1./9.*trace(b)*einsum('ij,kl',I,I) - \
            1./3.*(einsum('ij,kl',b,I) + einsum('ij,kl',I,b)) +\
            1./6.*trace(b)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)) )
        # VOLUMETRIC CHANGES
        if self.is_nearly_incompressible:
            H_Voigt += self.pressure*(einsum('ij,kl',I,I) - (einsum('ik,jl',I,I) + einsum('il,jk',I,I)))
        else:
            H_Voigt += self.kappa*((2.*J-1.)*einsum('ij,kl',I,I) - (J-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))

        H_Voigt = Voigt(H_Voigt,1)

        self.H_VoigtSize = H_Voigt.shape[0]

        return H_Voigt


    def CauchyStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        b = StrainTensors['b'][gcounter]

        mu = self.mu

        # ISOCHORIC PART
        stress = mu*J**(-5./3.)*(b - 1./3.*trace(b)*I)
        # VOLUMETRIC PART
        if self.is_nearly_incompressible:
            stress += self.pressure*I
        else:
            stress += self.kappa*(J-1.)*I

        return stress
