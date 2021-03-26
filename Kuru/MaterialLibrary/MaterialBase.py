from __future__ import print_function
import numpy as np
#from Florence.Utils import insensitive
from warnings import warn

# BASE CLASS FOR ALL MATERIAL MODELS - SHOULD NOT BE USED DIRECTLY
class Material(object):
    """Base class for all material models"""

    def __init__(self, mtype, ndim, energy_type="internal_energy", density=None,
        is_compressible=True, is_incompressible=False, is_nearly_incompressible=False,
        is_nonisotropic=True,is_anisotropic=False,is_transversely_isotropic=False,
        anisotropic_orientations=None, load_factor=None, **kwargs):


        # SAFETY CHECKS
        if not isinstance(mtype, str):
            raise TypeError("Type of material model should be given as a string")
        if not isinstance(energy_type, str):
            raise TypeError("Material energy can either be 'internal_energy' or 'enthalpy'")

        self.energy_type = energy_type

        # MATERIAL CONSTANTS
        self.rho = density

        # SET ALL THE OPTIONAL KEYWORD ARGUMENTS
        for i in kwargs.keys():
            if "__" not in i:
                setattr(self,i,kwargs[i])

        self.mtype = mtype
        self.ndim = ndim
        self.nvar = self.ndim

        self.H_Voigt = None

        if self.H_Voigt is not None:
            self.H_VoigtSize = self.H_Voigt.shape[0]

        # MIXED FORMULATION SELECTOR
        self.is_compressible = is_compressible
        self.is_nearly_incompressible = is_nearly_incompressible
        self.is_incompressible = is_incompressible

        self.is_anisotropic = is_anisotropic
        self.is_transversely_isotropic = is_transversely_isotropic
        self.is_nonisotropic = is_nonisotropic
        self.anisotropic_orientations = anisotropic_orientations

        # MIXED FORMULATION FOR VOLUME CHANGE
        self.pressure = 0.0

        self.has_low_level_dispatcher = False
        self.has_state_variables = False
        self.has_growth_remodeling = False

        # APPLY LOAD FACTOR TO DESPOSITION STRETCH
        self.load_factor = load_factor
        if self.load_factor is not None:
            self.factor_increment = self.load_factor[0]
        else:
            self.factor_increment = 1.0

    def ConnectivityOfMaterial(self,mesh):
        """
        Set conectivity between element_set and node_set for local nodes
        """
        elements = mesh.elements[self.element_set]
        aux = []
        for i in self.node_set:
            aux.append(np.where(elements==i))
        for j in range(self.node_set.shape[0]):
            elements[aux[j]]=j

        self.elements = np.array(elements, dtype=np.int64, copy=True)

    def GetNodeCommonality(self):
        """Finds the elements sharing a node.
            The return values are linked lists [list of numpy of arrays].
            Each numpy array within the list gives the elements that contain a given node.
            As a result the size of the linked list is nnode

            outputs:
                els:                        [list of numpy arrays] element numbers containing nodes
                pos:                        [list of numpy arrays] elemental positions of the nodes
                res_flat:                   [list of numpy arrays] position of nodes in the
                                            flattened element connectivity.
        """

#        self.__do_essential_memebers_exist__()

        elements = self.elements.ravel()
        idx_sort = np.argsort(elements)
        sorted_elements = elements[idx_sort]
        vals, idx_start = np.unique(sorted_elements, return_index=True)

        # Sets of indices
        flat_pos = np.split(idx_sort, idx_start[1:])
        els = np.split(idx_sort // int(self.elements.shape[1]), idx_start[1:])
        pos = np.split(idx_sort %  int(self.elements.shape[1]), idx_start[1:])

        return els, pos, flat_pos


    def MappingStateVariables(self,mesh,function_space,elem=0):
        """
        Function to map the field variables to the Gauss points in the element.
        """
        # Get index of the element into the element_set
        idx = np.where(self.element_set==elem)[0][0]

        # Get bases for the element
        Bases = function_space.Bases

        # Field variables at element nodes
        ElemStateVariables = self.state_variables[self.elements[idx,:],:]

        # Field variables at gauss points
        self.StateVariables = np.einsum('ij,ik->jk',Bases,ElemStateVariables)


