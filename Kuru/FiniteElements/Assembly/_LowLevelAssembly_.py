from warnings import warn
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

try:
    from ._LLADF__NeoHookean_2_ import _LLADF__NeoHookean_2_
    from ._LLADF__AnisotropicFungQuadratic_ import _LLADF__AnisotropicFungQuadratic_
    from ._LLADF__ArterialWallMixture_ import _LLADF__ArterialWallMixture_
    has_low_level_dispatcher = True
except ImportError:
    has_low_level_dispatcher = False
    warn("Cannot use low level dispatchers for Assembly")

__all__ = ['_LowLevelAssembly_']

def _LowLevelAssembly_(fem_solver, function_space, formulation, mesh, materials, Eulerx):

    prefix = "_LLADF__"

    # CHECK IF THE MATERIAL ASSEMBLY IS IMPLEMENTED
    for imat in range(len(materials)):
        assembly_func = prefix + type(materials[imat]).__name__ + "_"
        # CHECK IF LOW LEVEL ASSEMBLY EXISTS FOR MATERIAL
        ll_exists = False
        for key in globals().keys():
            if assembly_func == key:
                ll_exists = True
                break
        if ll_exists is False:
            raise NotImplementedError("Turning optimise option on for {} material is not supported yet. "
                "Consider 'optimise=False' for now".format(type(materials[imat]).__name__))

    nvar = formulation.nvar
    K = [[] for i in range(len(materials))]
    T = [[] for i in range(len(materials))]
    for imat in range(len(materials)):
        assembly_func = prefix + type(materials[imat]).__name__ + "_"
        # MAKE MESH DATA CONTIGUOUS
        mesh.ChangeType()
        # CALL LOW LEVEL ASSEMBLER
        if fem_solver.recompute_sparsity_pattern:
            I_stiffness, J_stiffness, V_stiffness, T[imat] = eval(assembly_func)(fem_solver,
                function_space, formulation, mesh, materials[imat], Eulerx)

            K[imat] = csr_matrix((V_stiffness,(I_stiffness,J_stiffness)),
                shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64)
        else:
            V_stiffness, T[imat] = eval(assembly_func)(fem_solver,
                function_space, formulation, mesh, materials[imat], Eulerx)

            K[imat] = csr_matrix((V_stiffness,fem_solver.indices,fem_solver.indptr),
                shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64)

    stiffness = K[0]
    TractionForce = T[0]
    if len(materials) > 1:
        for imat in range(len(materials)-1):
            stiffness += K[imat+1]
            TractionForce += T[imat+1]

    mass = []

    return stiffness, TractionForce, mass

