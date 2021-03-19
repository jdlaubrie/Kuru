import numpy as np

def KinematicMeasures(F,analysis_nature):
    """Computes the kinematic measures at all Gauss points for a given element.

    input:

        F:                  [ndarray (nofGauss x ndim x ndim)] Deformation gradient tensor evaluated at all integration points
        analysis_nature:       [str] type of analysis i.e. linear or non-linear 


    returns:

        StrainTensors:      [dict] a dictionary containing kinematic measures such as F, J, b etc evaluated at all integration points   
    """
    
    assert F.ndim == 3

    StrainTensors = {'F':F, 'J':np.linalg.det(F), 'b':np.einsum('ijk,ilk->ijl',F,F),
    'I':np.eye(F.shape[1],F.shape[1],dtype=np.float64)}

    # LINEARISED KINEMATICS
    # MATERIAL GRADIENT OF DISPLACEMENT
    StrainTensors['Gradu'] = F - StrainTensors['I']
    # SMALL STRAIN TENSOR IS THE LINEARISED VERSION OF GREEN-LAGRANGE STRAIN TENSOR
    StrainTensors['strain'] = 0.5*(StrainTensors['Gradu'] + np.einsum('ikj',StrainTensors['Gradu']))

        
    return StrainTensors

