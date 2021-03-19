from __future__ import print_function
import os, sys, gc
#from time import time
#from copy import deepcopy
import numpy as np
#import numpy.linalg as la
from warnings import warn
from Kuru import FunctionSpace, QuadratureRule
#from Florence.Base import JacobianError, IllConditionedError
#from Florence.Utils import PWD, RSWD

#from Florence.FunctionSpace import Tri
#from Florence.FunctionSpace import Tet
#from Florence.FunctionSpace import Quad, QuadES
#from Florence.FunctionSpace import Hex, HexES

from Kuru.FiniteElements.LocalAssembly.KinematicMeasures import *
from Kuru.FiniteElements.LocalAssembly._KinematicMeasures_ import _KinematicMeasures_
from Kuru import Mesh
#from Kuru import FEMSolver
#from Florence.MeshGeneration import vtk_writer
#from Florence.Utils import constant_camera_view

class PostProcess(object):
    """Post-process class for finite element solvers"""

    def __init__(self,ndim,nvar):

        self.domain_bases = None
        self.postdomain_bases = None
        self.boundary_bases = None
        self.ndim = ndim
        self.nvar = nvar
        self.analysis_type = None
        self.analysis_nature = None
        self.material_type = None

        self.is_scaledjacobian_computed = False
        self.is_material_anisotropic = False
        self.directions = None

        self.mesh = None
        self.sol = None
        self.recovered_fields = None

        self.formulation = None
        self.material = None
        self.fem_solver = None

        self.parallel_model = None
        self.ncpu = None

    def SetBases(self,domain=None,postdomain=None,boundary=None):
        """Sets bases for all integration points for 'domain', 'postdomain' or 'boundary'
        """

        if domain is None and postdomain is None and boundary is None:
            warn("Nothing to be set")

        self.domain_bases = domain
        self.postdomain_bases = postdomain
        self.boundary_bases = boundary

    def SetAnalysis(self,analysis_type=None, analysis_nature=None):
        self.analysis_type = analysis_type
        self.analysis_nature = analysis_nature

    def SetMesh(self,mesh):
        """Set initial (undeformed) mesh"""
        self.mesh = mesh

    def SetSolution(self,sol):
        self.sol = sol

    def SetFormulation(self,formulation):
        self.formulation = formulation

    def SetMaterial(self,materials):
        self.materials = materials

    def SetFEMSolver(self,fem_solver):
        self.fem_solver = fem_solver

    def SetGrowthRemodeling(self,gr_variables):
        self.gr_variables = gr_variables

    def NodeStressRecovery(self, mynode=0, imat=0, steps=None):

        """
            steps:          [list,np.1darray] for which time steps/increments the data should
                            be recovered
        """

        if self.mesh is None:
            raise ValueError("Mesh not set for post-processing")
        if self.sol is None:
            raise ValueError("Solution not set for post-processing")
        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.materials is None:
            raise ValueError("materials not set for post-processing")
        if self.fem_solver is None:
            raise ValueError("FEM solver not set for post-processing")

        if self.sol.shape[1] > self.nvar:
            return

        det = np.linalg.det
        inv = np.linalg.inv

        mesh = self.mesh
        fem_solver = self.fem_solver
        formulation = self.formulation
        material = self.materials[imat]

        if not mynode in material.node_set:
            raise ValueError("Node {} is not in material {}".format(mynode,imat))

        # GET THE UNDERLYING LINEAR MESH
        # lmesh = mesh.GetLinearMesh()
        C = mesh.InferPolynomialDegree() - 1
        ndim = mesh.InferSpatialDimension()

        elements = mesh.elements
        points = mesh.points
        nelem = elements.shape[0]; npoint = points.shape[0]
        nodeperelem = elements.shape[1]

        # GET QUADRATURE
        norder = 2*C
        if norder == 0:
            norder=1
        # quadrature = QuadratureRule(qtype="gauss", norder=norder, mesh_type=mesh.element_type, optimal=3)
        # Domain = FunctionSpace(mesh, quadrature, p=C+1)
        Domain = FunctionSpace(mesh, p=C+1, evaluate_at_nodes=True)
        Jm = Domain.Jm
        AllGauss = Domain.AllGauss
        Bases = Domain.Bases

        # requires_geometry_update = fem_solver.requires_geometry_update
        requires_geometry_update = True # ALWAYS TRUE FOR THIS ROUTINE
        TotalDisp = self.sol[:,:]

        if hasattr(self,"number_of_time_increments"):
            TimeIncrement = self.number_of_time_increments
        else:
            TimeIncrement = fem_solver.number_of_time_increments
        increments = range(TimeIncrement)
        if steps is not None:
            TimeIncrement = len(steps)
            increments = steps

        material_node = np.where(material.node_set==mynode)[0][0]
        # COMPUTE THE COMMON/NEIGHBOUR NODES ONCE
        Elss, Poss = mesh.GetNodeCommonality()[:2]
        Pos = Poss[mynode]

        MainDict = {}
        MainDict['F'] = np.zeros((TimeIncrement,ndim,ndim))
        MainDict['CauchyStress'] = np.zeros((TimeIncrement,ndim,ndim))
        if self.gr_variables is not None:
            MainDict['FibreStress'] = np.zeros((TimeIncrement,5))

        material.ConnectivityOfMaterial(mesh)
        Elsm, Posm = material.GetNodeCommonality()[:2]
        Elm = Elsm[material_node]
        ncommon_nodes_m = Elm.shape[0]

        for incr, Increment in enumerate(increments):
            if TotalDisp.ndim == 3:
                Eulerx = points + TotalDisp[:,:ndim,Increment]
            else:
                Eulerx = points + TotalDisp[:,:ndim]

            F = np.zeros((ncommon_nodes_m,nodeperelem,ndim,ndim))
            CauchyStress = np.zeros((ncommon_nodes_m,ndim,ndim))
            if self.gr_variables is not None:
                FibreStress = np.zeros((ncommon_nodes_m,5))

            # LOOP OVER ELEMENTS
            for i in range(ncommon_nodes_m):
                ielem = Elm[i]
                elem = material.element_set[ielem]
                # GET THE FIELDS AT THE ELEMENT LEVEL
                LagrangeElemCoords = points[elements[elem,:],:]
                EulerELemCoords = Eulerx[elements[elem,:],:]
                # GROWTH-REMODELING VALUES FOR THIS ELEMENT
                if material.has_state_variables:
                    material.state_variables[:,9:21] = self.gr_variables[imat][:,:,Increment]
                    material.MappingStateVariables(mesh,Domain,elem)

                if material.has_low_level_dispatcher:

                    CauchyStressAux = np.zeros((nodeperelem,ndim,ndim))
                    FibreStressAux = np.zeros((nodeperelem,5))
                    # GET LOCAL KINEMATICS
                    SpatialGradient, F[i,:,:,:], detJ, dV = _KinematicMeasures_(Jm, AllGauss[:,0],
                            LagrangeElemCoords, EulerELemCoords, requires_geometry_update)
                    # PARAMETERS FOR INCOMPRESSIBILITY (MEAN DILATATION METHOD HU-WASHIZU)
                    if material.is_incompressible:
                        MaterialVolume = np.sum(dV)
                        if material.has_growth_remodeling:
                            dve = np.true_divide(detJ,material.StateVariables[:,20])
                            CurrentVolume = np.sum(dve)
                        else:
                            CurrentVolume = np.sum(detJ)
                        material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                    # COMPUTE WORK-CONJUGATES AND HESSIAN AT THIS GAUSS POINT
                    counter = Pos[i]
                    CauchyStressAux[:,:], _ = material.KineticMeasures(F[i,:,:,:],elem=elem)
                    CauchyStress[i,:] = CauchyStressAux[counter,:]
                    if self.gr_variables is not None:
                        FibreStressAux[:,:], _ = material.LLConstituentStress(F[i,:,:,:],elem=elem)
                        FibreStress[i,:] = FibreStressAux[counter,:]

                else:
                    # GAUSS LOOP IN VECTORISED FORM
                    ParentGradientX = np.einsum('ijk,jl->kil', Jm, LagrangeElemCoords)
                    # MATERIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla_0 (N)]
                    MaterialGradient = np.einsum('ijk,kli->ijl', inv(ParentGradientX), Jm)
                    # DEFORMATION GRADIENT TENSOR [\vec{x} \otimes \nabla_0 (N)]
                    F[i,:,:,:] = np.einsum('ij,kli->kjl', EulerELemCoords, MaterialGradient)
                    # COMPUTE REMAINING KINEMATIC MEASURES
                    StrainTensors = KinematicMeasures(F[i,:,:,:], fem_solver.analysis_nature)

                    # GEOMETRY UPDATE IS A MUST
                    # MAPPING TENSOR [\partial\vec{X}/ \partial\vec{\varepsilon} (ndim x ndim)]
                    ParentGradientx = np.einsum('ijk,jl->kil',Jm,EulerELemCoords)
                    # SPATIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla (N)]
                    SpatialGradient = np.einsum('ijk,kli->ilj',inv(ParentGradientx),Jm)
                    # COMPUTE ONCE detJ (GOOD SPEEDUP COMPARED TO COMPUTING TWICE)
                    detJ = np.einsum('i,i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)),
                        np.abs(StrainTensors['J']))

                    # COMPUTE PARAMETERS FOR MEAN DILATATION METHOD, IT NEEDS TO BE BEFORE COMPUTE HESSIAN AND STRESS
                    if material.is_incompressible:
                        dV = np.einsum('i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)))
                        MaterialVolume = np.sum(dV)
                        if material.has_growth_remodeling:
                            dve = np.true_divide(detJ,material.StateVariables[:,20])
                            CurrentVolume = np.sum(dve)
                        else:
                            CurrentVolume = np.sum(detJ)
                        material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                    counter = Pos[i]
                    CauchyStress[i,:] = material.CauchyStress(StrainTensors,elem,counter)
                    if self.gr_variables is not None:
                        FibreStress[i,:],_ = material.ConstituentStress(StrainTensors,elem,counter)

            for i in range(ncommon_nodes_m):
                MainDict['F'][incr,:,:] += F[i,Pos[i],:,:]
                MainDict['CauchyStress'][incr,:,:] += CauchyStress[i,:,:]
                if self.gr_variables is not None:
                    MainDict['FibreStress'][incr,:] += FibreStress[i,:]
            # AVERAGE OUT
            MainDict['F'][incr,:,:] /= ncommon_nodes_m
            MainDict['CauchyStress'][incr,:,:] /= ncommon_nodes_m
            if self.gr_variables is not None:
                MainDict['FibreStress'][incr,:] /= ncommon_nodes_m


        self.node_recovered_fields = MainDict
        return

    def StressRecovery(self, steps=None, average_derived_quantities=True, time_problem=True):

        """
            steps:          [list,np.1darray] for which time steps/increments the data should
                            be recovered
        """

        if self.mesh is None:
            raise ValueError("Mesh not set for post-processing")
        if self.sol is None:
            raise ValueError("Solution not set for post-processing")
        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.materials is None:
            raise ValueError("materials not set for post-processing")
        if self.fem_solver is None:
            raise ValueError("FEM solver not set for post-processing")

        if self.sol.shape[1] > self.nvar:
            return

        det = np.linalg.det
        inv = np.linalg.inv

        mesh = self.mesh
        fem_solver = self.fem_solver
        formulation = self.formulation
        materials = self.materials

        # GET THE UNDERLYING LINEAR MESH
        # lmesh = mesh.GetLinearMesh()
        C = mesh.InferPolynomialDegree() - 1
        ndim = mesh.InferSpatialDimension()

        elements = mesh.elements
        points = mesh.points
        nelem = elements.shape[0]; npoint = points.shape[0]
        nodeperelem = elements.shape[1]

        # GET QUADRATURE
        norder = 2*C
        if norder == 0:
            norder=1
        # quadrature = QuadratureRule(qtype="gauss", norder=norder, mesh_type=mesh.element_type, optimal=3)
        # Domain = FunctionSpace(mesh, quadrature, p=C+1)
        Domain = FunctionSpace(mesh, p=C+1, evaluate_at_nodes=True)
        Jm = Domain.Jm
        AllGauss = Domain.AllGauss
        Bases = Domain.Bases

        # requires_geometry_update = fem_solver.requires_geometry_update
        requires_geometry_update = True # ALWAYS TRUE FOR THIS ROUTINE
        TotalDisp = self.sol[:,:]

        if time_problem is True:
            if hasattr(self,"number_of_time_increments"):
                TimeIncrement = self.number_of_time_increments
            else:
                TimeIncrement = fem_solver.number_of_time_increments
            increments = range(TimeIncrement)
            if steps is not None:
                TimeIncrement = len(steps)
                increments = steps
            MainDict = {}
            MainDict['F'] = np.zeros((TimeIncrement,npoint,ndim,ndim))
            MainDict['CauchyStress'] = [[] for i in range(len(materials))]
            MainDict['FibreStress'] = [[] for i in range(len(materials))]
            for imat in range(len(materials)):
                materials[imat].ConnectivityOfMaterial(mesh)
                MainDict['CauchyStress'][imat] = np.zeros((TimeIncrement,materials[imat].node_set.shape[0],ndim,ndim))
                if self.gr_variables is not None:
                    MainDict['FibreStress'][imat] = np.zeros((TimeIncrement,materials[imat].node_set.shape[0],5))
        else:
            if hasattr(self,"number_of_load_increments"):
                LoadIncrement = self.number_of_time_increments
            else:
                LoadIncrement = fem_solver.number_of_load_increments
            increments = range(LoadIncrement)
            if steps is not None:
                LoadIncrement = len(steps)
                increments = steps
            MainDict = {}
            MainDict['F'] = np.zeros((LoadIncrement,npoint,ndim,ndim))
            MainDict['CauchyStress'] = [[] for i in range(len(materials))]
            MainDict['FibreStress'] = [[] for i in range(len(materials))]
            for imat in range(len(materials)):
                materials[imat].ConnectivityOfMaterial(mesh)
                MainDict['CauchyStress'][imat] = np.zeros((LoadIncrement,materials[imat].node_set.shape[0],ndim,ndim))
                if self.gr_variables is not None:
                    MainDict['FibreStress'][imat] = np.zeros((LoadIncrement,materials[imat].node_set.shape[0],5))

        # COMPUTE THE COMMON/NEIGHBOUR NODES ONCE
        all_nodes = np.unique(elements)
        Elss, Poss = mesh.GetNodeCommonality()[:2]

        for incr, Increment in enumerate(increments):

            if TotalDisp.ndim == 3:
                Eulerx = points + TotalDisp[:,:ndim,Increment]
            else:
                Eulerx = points + TotalDisp[:,:ndim]

            F = np.zeros((nelem,nodeperelem,ndim,ndim))

            for imat in range(len(materials)):
                material = materials[imat]
                Elsm, Posm = material.GetNodeCommonality()[:2]
                CauchyStress = np.zeros((material.element_set.shape[0],nodeperelem,ndim,ndim))
                if self.gr_variables is not None:
                    FibreStress = np.zeros((material.element_set.shape[0],nodeperelem,5))
                # LOOP OVER ELEMENTS
                for ielem in range(material.element_set.shape[0]):
                    elem = material.element_set[ielem]
                    # GET THE FIELDS AT THE ELEMENT LEVEL
                    LagrangeElemCoords = points[elements[elem,:],:]
                    EulerELemCoords = Eulerx[elements[elem,:],:]
                    # GROWTH-REMODELING VALUES FOR THIS ELEMENT
                    if material.has_state_variables:
                        if self.gr_variables is None:
                            material.MappingStateVariables(mesh,Domain,elem)
                        elif self.gr_variables is not None:
                            material.state_variables[:,9:21] = self.gr_variables[imat][:,:,Increment]
                            material.MappingStateVariables(mesh,Domain,elem)

                    if material.has_low_level_dispatcher:

                        # GET LOCAL KINEMATICS
                        SpatialGradient, F[elem,:,:,:], detJ, dV = _KinematicMeasures_(Jm, AllGauss[:,0],
                                LagrangeElemCoords, EulerELemCoords, requires_geometry_update)
                        # PARAMETERS FOR INCOMPRESSIBILITY (MEAN DILATATION METHOD HU-WASHIZU)
                        if material.is_incompressible:
                            MaterialVolume = np.sum(dV)
                            if material.has_growth_remodeling:
                                dve = np.true_divide(detJ,material.StateVariables[:,20])
                                CurrentVolume = np.sum(dve)
                            else:
                                CurrentVolume = np.sum(detJ)
                            material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                        # COMPUTE WORK-CONJUGATES AND HESSIAN AT THIS GAUSS POINT
                        CauchyStress[ielem,:,:], _ = material.KineticMeasures(F[elem,:,:,:],elem=elem)
                        if self.gr_variables is not None and material.has_state_variables:
                            FibreStress[ielem,:,:], _ = material._ConstituentMeasures_(F[elem,:,:,:],elem=elem)

                    else:
                        # GAUSS LOOP IN VECTORISED FORM
                        ParentGradientX = np.einsum('ijk,jl->kil', Jm, LagrangeElemCoords)
                        # MATERIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla_0 (N)]
                        MaterialGradient = np.einsum('ijk,kli->ijl', inv(ParentGradientX), Jm)
                        # DEFORMATION GRADIENT TENSOR [\vec{x} \otimes \nabla_0 (N)]
                        F[elem,:,:,:] = np.einsum('ij,kli->kjl', EulerELemCoords, MaterialGradient)
                        # COMPUTE REMAINING KINEMATIC MEASURES
                        StrainTensors = KinematicMeasures(F[elem,:,:,:], fem_solver.analysis_nature)

                        # GEOMETRY UPDATE IS A MUST
                        # MAPPING TENSOR [\partial\vec{X}/ \partial\vec{\varepsilon} (ndim x ndim)]
                        ParentGradientx = np.einsum('ijk,jl->kil',Jm,EulerELemCoords)
                        # SPATIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla (N)]
                        SpatialGradient = np.einsum('ijk,kli->ilj',inv(ParentGradientx),Jm)
                        # COMPUTE ONCE detJ (GOOD SPEEDUP COMPARED TO COMPUTING TWICE)
                        detJ = np.einsum('i,i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)),
                            np.abs(StrainTensors['J']))

                        # COMPUTE PARAMETERS FOR MEAN DILATATION METHOD, IT NEEDS TO BE BEFORE COMPUTE HESSIAN AND STRESS
                        if material.is_incompressible:
                            dV = np.einsum('i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)))
                            MaterialVolume = np.sum(dV)
                            if material.has_growth_remodeling:
                                dve = np.true_divide(detJ,material.StateVariables[:,20])
                                CurrentVolume = np.sum(dve)
                            else:
                                CurrentVolume = np.sum(detJ)
                            material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                        # LOOP OVER GAUSS POINTS
                        for counter in range(AllGauss.shape[0]):
                            CauchyStress[ielem,counter,:] = material.CauchyStress(StrainTensors,elem,counter)
                            if self.gr_variables is not None and material.has_state_variables:
                                FibreStress[ielem,counter,:],_ = material.ConstituentMeasures(StrainTensors,elem,counter)

                if average_derived_quantities:
                    for inode in range(material.node_set.shape[0]):
                        Els, Pos = Elsm[inode], Posm[inode]
                        ncommon_nodes = Els.shape[0]
                        for uelem in range(ncommon_nodes):
                            MainDict['CauchyStress'][imat][incr,inode,:,:] += CauchyStress[Els[uelem],Pos[uelem],:,:]
                            if self.gr_variables is not None:
                                MainDict['FibreStress'][imat][incr,inode,:] += FibreStress[Els[uelem],Pos[uelem],:]
                        # AVERAGE OUT
                        MainDict['CauchyStress'][imat][incr,inode,:,:] /= ncommon_nodes
                        if self.gr_variables is not None:
                            MainDict['FibreStress'][imat][incr,inode,:] /= ncommon_nodes
                else:
                    for inode in range(material.node_set.shape[0]):
                        Els, Pos = Elsm[inode], Posm[inode]
                        ncommon_nodes = Els.shape[0]
                        uelem = 0
                        MainDict['CauchyStress'][imat][incr,inode,:,:] = CauchyStress[Els[uelem],Pos[uelem],:,:]
                        if self.gr_variables is not None:
                            MainDict['FibreStress'][imat][incr,inode,:] = FibreStress[Els[uelem],Pos[uelem],:]

            if average_derived_quantities:
                for inode in all_nodes:
                    Els, Pos = Elss[inode], Poss[inode]
                    ncommon_nodes = Els.shape[0]
                    for uelem in range(ncommon_nodes):
                        MainDict['F'][incr,inode,:,:] += F[Els[uelem],Pos[uelem],:,:]
                    # AVERAGE OUT
                    MainDict['F'][incr,inode,:,:] /= ncommon_nodes
            else:
                for inode in all_nodes:
                    Els, Pos = Elss[inode], Poss[inode]
                    ncommon_nodes = Els.shape[0]
                    uelem = 0
                    MainDict['F'][incr,inode,:,:] = F[Els[uelem],Pos[uelem],:,:]


        self.recovered_fields = MainDict
        return

    def AverageDeformationGradient(self, element_sets, fibre_direction):

        """
            steps:          [list,np.1darray] for which time steps/increments the data should
                            be recovered
        """

        if self.mesh is None:
            raise ValueError("Mesh not set for post-processing")
        if self.sol is None:
            raise ValueError("Solution not set for post-processing")
        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.materials is None:
            raise ValueError("materials not set for post-processing")
        if self.fem_solver is None:
            raise ValueError("FEM solver not set for post-processing")

        if self.sol.shape[1] > self.nvar:
            return

        det = np.linalg.det
        inv = np.linalg.inv

        mesh = self.mesh
        fem_solver = self.fem_solver
        formulation = self.formulation
        materials = self.materials

        # GET THE UNDERLYING LINEAR MESH
        # lmesh = mesh.GetLinearMesh()
        C = mesh.InferPolynomialDegree() - 1
        ndim = mesh.InferSpatialDimension()

        elements = mesh.elements
        points = mesh.points
        nelem = elements.shape[0]; npoint = points.shape[0]
        nodeperelem = elements.shape[1]

        # GET QUADRATURE
        norder = 2*C
        if norder == 0:
            norder=1
        # quadrature = QuadratureRule(qtype="gauss", norder=norder, mesh_type=mesh.element_type, optimal=3)
        # Domain = FunctionSpace(mesh, quadrature, p=C+1)
        Domain = FunctionSpace(mesh, p=C+1, evaluate_at_nodes=True)
        Jm = Domain.Jm
        AllGauss = Domain.AllGauss
        Bases = Domain.Bases

        # requires_geometry_update = fem_solver.requires_geometry_update
        requires_geometry_update = True # ALWAYS TRUE FOR THIS ROUTINE
        TotalDisp = self.sol[:,:]


        # COMPUTE THE COMMON/NEIGHBOUR NODES ONCE
        all_nodes = np.unique(elements)
        Elss, Poss = mesh.GetNodeCommonality()[:2]

        I = np.eye(3,3,dtype=np.float64)

        if TotalDisp.ndim == 3:
            Eulerx = points + TotalDisp[:,:ndim,-1]
        else:
            Eulerx = points + TotalDisp[:,:ndim]


        DeformationGradient = np.zeros((npoint,ndim,ndim))
        F_local = np.zeros((nelem,ndim,ndim))

        # LOOP OVER ELEMENTS
        for elem in range(nelem):
            # find the set of the element
            iset = -1
            for imat in range(len(materials)):
                if elem in element_sets[imat]:
                    iset = imat
                    break
            if iset is -1:
                raise ValueError("Set index is not well recognized")

            material = materials[iset]
            # GET THE FIELDS AT THE ELEMENT LEVEL
            LagrangeElemCoords = points[elements[elem,:],:]
            EulerELemCoords = Eulerx[elements[elem,:],:]

            if material.has_low_level_dispatcher:
                # GET LOCAL KINEMATICS
                SpatialGradient, F, detJ, dV = _KinematicMeasures_(Jm, AllGauss[:,0],
                        LagrangeElemCoords, EulerELemCoords, requires_geometry_update)

            else:
                # GAUSS LOOP IN VECTORISED FORM
                ParentGradientX = np.einsum('ijk,jl->kil', Jm, LagrangeElemCoords)
                # MATERIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla_0 (N)]
                MaterialGradient = np.einsum('ijk,kli->ijl', inv(ParentGradientX), Jm)
                # DEFORMATION GRADIENT TENSOR [\vec{x} \otimes \nabla_0 (N)]
                F = np.einsum('ij,kli->kjl', EulerELemCoords, MaterialGradient)

            # Directional tensor for element
            Normal = fibre_direction[elem][0][:,None]
            Normal = np.dot(I,Normal)[:,0]
            Tangential = fibre_direction[elem][1][:,None]
            Tangential = np.dot(I,Tangential)[:,0]
            Axial = fibre_direction[elem][2][:,None]
            Axial = np.dot(I,Axial)[:,0]
            Rotation = np.eye(3,3,dtype=np.float64)
            for i in range(3):
                Rotation[0,i] = Normal[i]
                Rotation[1,i] = Tangential[i]
                Rotation[2,i] = Axial[i]

            for counter in range(AllGauss.shape[0]):
                # cylindric world
                F_local[elem,:,:] += np.dot(Rotation,np.dot(F[counter,:,:],Rotation.T))
            F_local[elem,:,:] /= AllGauss.shape[0]

        for inode in all_nodes:
            Els, Pos = Elss[inode], Poss[inode]
            ncommon_nodes = Els.shape[0]
            for uelem in range(ncommon_nodes):
                DeformationGradient[inode,:,:] += F_local[Els[uelem],:,:]
            # AVERAGE OUT
            DeformationGradient[inode,:,:] /= ncommon_nodes

        self.average_deformation = DeformationGradient
        return


    def QuantityNamer(self, num, print_name=True):
        """Returns the quantity (for augmented solution i.e. primary and recovered variables)
            name given its number (from numbering order)

            Variables (quantity) numbering:

                quantity    mechanics 2D    mechanics 3D    growth_remodeling 2D    growth_remodeling 3D
                ----------------------------------------------------------------------------------------
                0           ux              ux              ux                      ux
                ----------------------------------------------------------------------------------------
                1           uy              uy              uy                      uy
                ----------------------------------------------------------------------------------------
                2           F_xx            uz              F_xx                    uz
                ----------------------------------------------------------------------------------------
                3           F_xy            F_xx            F_xy                    F_xx
                ----------------------------------------------------------------------------------------
                4           F_yx            F_xy            F_yx                    F_xy
                ----------------------------------------------------------------------------------------
                5           F_yy            F_xz            F_yy                    F_xz
                ----------------------------------------------------------------------------------------
                6           H_xx            F_yx            J                       F_yx
                ----------------------------------------------------------------------------------------
                7           H_xy            F_yy            C_xx                    F_yy
                ----------------------------------------------------------------------------------------
                8           H_yx            F_yz            C_xy                    F_yz
                ----------------------------------------------------------------------------------------
                9           H_yy            F_zx            C_yy                    F_zx
                ----------------------------------------------------------------------------------------
                10          J               F_zy            detC                    F_zy
                ----------------------------------------------------------------------------------------
                11          C_xx            F_zz            e_xx                    F_zz
                ----------------------------------------------------------------------------------------
                12          C_xy            H_xx            e_xy                    J
                ----------------------------------------------------------------------------------------
                13          C_yy            H_xy            e_yy                    C_xx
                ----------------------------------------------------------------------------------------
                14          G_xx            H_xz            e_hyd                   C_xy
                ----------------------------------------------------------------------------------------
                15          G_xy            H_yx            e_VM                    C_xz
                ----------------------------------------------------------------------------------------
                16          G_yy            H_yy            s_xx                    C_yy
                ----------------------------------------------------------------------------------------
                17          detC            H_yz            s_xy                    C_yz
                ----------------------------------------------------------------------------------------
                18          S_xx            H_zx            s_yy                    C_zz
                ----------------------------------------------------------------------------------------
                19          S_xy            H_zy            p_hyd                   detC
                ----------------------------------------------------------------------------------------
                20          S_yy            H_zz            s_VM                    e_xx
                ----------------------------------------------------------------------------------------
                21          p_hyd           J               s_m                     e_xy
                ----------------------------------------------------------------------------------------
                22          None            C_xx            s_c1                    e_xz
                ----------------------------------------------------------------------------------------
                23          None            C_xy            s_c2                    e_yy
                ----------------------------------------------------------------------------------------
                24          None            C_xz            s_c3                    e_yz
                ----------------------------------------------------------------------------------------
                25          None            C_yy            s_c4                    e_zz
                ----------------------------------------------------------------------------------------
                26          None            C_yz            den_e                   e_hyd
                ----------------------------------------------------------------------------------------
                27          None            C_zz            den_m                   e_VM
                ----------------------------------------------------------------------------------------
                28          None            G_xx            den_c                   s_xx
                ----------------------------------------------------------------------------------------
                29          None            G_xy            growth                  s_xy
                ----------------------------------------------------------------------------------------
                30          None            G_xz            rem_m                   s_xz
                ----------------------------------------------------------------------------------------
                31          None            G_yy            rem_c1                  s_yy
                ----------------------------------------------------------------------------------------
                32          None            G_yz            rem_c2                  s_yz
                ----------------------------------------------------------------------------------------
                33          None            G_zz            rem_c3                  s_zz
                ----------------------------------------------------------------------------------------
                34          None            detC            rem_c4                  p_hyd
                ----------------------------------------------------------------------------------------
                35          None            S_xx            None                    s_VM
                ----------------------------------------------------------------------------------------
                36          None            S_xy            None                    s_m
                ----------------------------------------------------------------------------------------
                37          None            S_xz            None                    s_c1
                ----------------------------------------------------------------------------------------
                38          None            S_yy            None                    s_c2
                ----------------------------------------------------------------------------------------
                39          None            S_yz            None                    s_c3
                ----------------------------------------------------------------------------------------
                40          None            S_zz            None                    s_c4
                ----------------------------------------------------------------------------------------
                41          None            p_hyd           None                    den_e
                ----------------------------------------------------------------------------------------
                42          None            s_VM            None                    den_m
                ----------------------------------------------------------------------------------------
                43          None            None            None                    den_c
                ----------------------------------------------------------------------------------------
                44          None            None            None                    growth
                ----------------------------------------------------------------------------------------
                45          None            None            None                    rem_m
                ----------------------------------------------------------------------------------------
                46          None            None            None                    rem_c1
                ----------------------------------------------------------------------------------------
                47          None            None            None                    rem_c2
                ----------------------------------------------------------------------------------------
                48          None            None            None                    rem_c3
                ----------------------------------------------------------------------------------------
                49          None            None            None                    rem_c4
                ----------------------------------------------------------------------------------------





            where S represents Cauchy stress tensor, E the strain tenosr, {S_m,S_ci} the fibre stresses,
            {den_e,den_m,den_ci} the densities and {rem_m,rem_ci} the fibre remodelig
        """
        namer = None
        if num > 49:
            if print_name:
                print('Quantity corresponds to ' + str(namer))
            return namer

        lines = []
        with open(__file__) as f:
            lines.append(f.readlines())
        lines = lines[0]

        line_number = len(lines)+1

        for counter, line in enumerate(lines):
            line = line.strip()
            if "quantity" in line and "mechanics" in line and "2D" in line and "3D" in line:
                line_number = counter
            if counter > line_number+1 and counter < line_number+105:
                spl = list(filter(None, line.split(" ")))
                if spl[0] == str(num):
                    if self.nvar == 2 and self.ndim==2 and self.gr_variables is None:
                        namer = spl[-4]
                    elif self.nvar == 2 and self.ndim==2 and self.gr_variables is not None:
                        namer = spl[-2]
                    elif self.nvar == 3 and self.ndim==3 and self.gr_variables is None:
                        namer = spl[-3]
                    elif self.nvar == 3 and self.ndim==3 and self.gr_variables is not None:
                        namer = spl[-1]
                    break


        if print_name:
            print('Quantity corresponds to ' + str(namer))

        if "ux" in namer:
            namer = "u_x"
        elif "uy" in namer:
            namer = "u_y"
        elif "uz" in namer:
            namer = "u_z"
        elif "phi" in namer:
            namer = "\phi"
        return namer

    def GetAugmentedSolution(self, steps=None, average_derived_quantities=True, parallelise=False, time_problem=True):
        """This function modifies self.sol to augmented_sol and returns the augmented solution
            augmented_sol
        """

        if self.sol.shape[1] > self.nvar:
            return self.sol

        import inspect
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)[1][3]
        if calframe != "ParallelGetAugmentedSolution":
            print("Computing recovered quantities. This is going to take some time...")

        if self.fem_solver is None:
            raise ValueError("FEM solver not set for post-processing")
        if steps is None:
            if self.fem_solver.number_of_time_increments == 1:
                parallelise = False
        else:
            if len(steps) == 1:
                parallelise = False

        if parallelise:
            from multiprocessing import Pool, cpu_count
            from contextlib import closing
            if self.ncpu is None:
                self.ncpu = cpu_count()
            if self.parallel_model is None:
                self.parallel_model = "pool"

            if steps is not None:
                increments = steps
            else:
                increments = range(self.fem_solver.number_of_time_increments)
            increments = np.array(increments).flatten()
            partitioned_steps = np.array_split(increments,self.ncpu)

            pps = []
            for ip in range(len(partitioned_steps)):
                pp = deepcopy(self)
                pp.number_of_time_increments = len(partitioned_steps[ip])
                zipper_object = ParallelGetAugmentedSolutionZipper(pp,
                    steps=partitioned_steps[ip], average_derived_quantities=average_derived_quantities)
                pps.append(zipper_object)

            # Pool
            if self.parallel_model == "pool":
                with closing(Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelGetAugmentedSolution,pps)
                    pool.terminate()

            # Thread Pool
            elif self.parallel_model == "thread_pool":
                import multiprocessing.dummy
                with closing(multiprocessing.dummy.Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelGetAugmentedSolution,pps)
                    pool.terminate()

            else:
                raise ValueError("Parallel model not understood")

            # Serial
            # res = []
            # for ip in range(len(partitioned_steps)):
            #     res.append(ParallelGetAugmentedSolution(pps[ip]))


            ndim = self.formulation.ndim
            fields = self.formulation.fields
            nnode = self.mesh.points.shape[0]
            if fields == "mechanics" and ndim == 2:
                augmented_sol = np.zeros((nnode,22,len(increments)),dtype=np.float64)
            elif fields == "mechanics" and ndim == 3:
                augmented_sol = np.zeros((nnode,42,len(increments)),dtype=np.float64)

            for ip in range(len(partitioned_steps)):
                augmented_sol[:,:,partitioned_steps[ip]] = res[ip]

            self.sol = augmented_sol
            return augmented_sol


        # GET RECOVERED VARIABLES ALL VARIABLE CHECKS ARE DONE IN STRESS RECOVERY
        self.StressRecovery(steps=steps, average_derived_quantities=average_derived_quantities, time_problem=time_problem)

        ndim = self.formulation.ndim
        fields = self.formulation.fields
        nnode = self.mesh.points.shape[0]
        if hasattr(self, "number_of_time_increments"):
            increments = self.number_of_time_increments
        else:
            if self.sol.ndim == 3:
                increments = self.sol.shape[2]
            else:
                increments = 1
            if steps is not None:
                increments = len(steps)
            else:
                if increments != self.fem_solver.number_of_time_increments and increments != self.fem_solver.number_of_load_increments:
                    raise ValueError("Incosistent number of time/load increments between FEMSolver and provided solution")


        # Get tensors
        if ndim == 2:
            I=np.eye(2)
        elif ndim == 3:
            I=np.eye(3)

        nnodes = 0
        for imat in range(len(self.materials)):
            nnodes += self.materials[imat].node_set.shape[0]

        if fields == "mechanics" and ndim == 2 and self.gr_variables is None:
            augmented_sol = np.zeros((nnodes,22,increments),dtype=np.float64)
        elif fields == "mechanics" and ndim == 3 and self.gr_variables is None:
            augmented_sol = np.zeros((nnodes,43,increments),dtype=np.float64)
        elif fields == "mechanics" and ndim == 2 and self.gr_variables is not None:
            augmented_sol = np.zeros((nnodes,35,increments),dtype=np.float64)
        elif fields == "mechanics" and ndim == 3 and self.gr_variables is not None:
            augmented_sol = np.zeros((nnodes,50,increments),dtype=np.float64)

        nstart = 0
        nend = 0
        for imat in range(len(self.materials)):
            nnode_set = self.materials[imat].node_set.shape[0]
            nend += nnode_set
            F = self.recovered_fields['F'][:,self.materials[imat].node_set,:,:]
            J = np.linalg.det(F)
            C = np.einsum('ijlk,ijkm->ijlm',np.einsum('ijlk',F),F)
            detC = J**2
            Cauchy = self.recovered_fields['CauchyStress'][imat]
            p_hyd = 1./3.*np.einsum('ijkk',Cauchy)
            # Get tensor for just mechanics fields or for growth and remodeling
            if self.gr_variables is None:
                H = np.einsum('ij,ijlk->ijkl',J,np.linalg.inv(F))
                G = np.einsum('ijlk,ijkm->ijlm',np.einsum('ijlk',H),H)
                s_d = Cauchy - np.einsum('ij,kl->ijkl',p_hyd,I)
                s_VM = np.sqrt(3./2.*np.einsum('ijkl,ijkl->ij',s_d,s_d))
                # Reshape
                H = np.einsum('lijk',H).reshape(nnode_set,ndim**2,increments)
                G = np.einsum('lijk',G).reshape(nnode_set,ndim**2,increments)
                s_VM = np.einsum('ji',s_VM).reshape(nnode_set,increments)
            else:
                b_inv = np.einsum('ijlk,ijkm->ijlm',np.linalg.inv(F),np.linalg.inv(F))
                e = 1./2.*(I-b_inv)
                e_hyd = 1./3.*np.einsum('ijkk',e)
                e_d = e - np.einsum('ij,kl->ijkl',e_hyd,I)
                e_VM = np.sqrt(2./3.*np.einsum('ijkl,ijkl->ij',e_d,e_d))
                s_d = Cauchy - np.einsum('ij,kl->ijkl',p_hyd,I)
                s_VM = np.sqrt(3./2.*np.einsum('ijkl,ijkl->ij',s_d,s_d))
                FibreStress = self.recovered_fields['FibreStress'][imat]
                collagen_density = np.einsum('ijk->ik',self.gr_variables[imat][:,7:11,:])
                # Reshape
                e = np.einsum('lijk',e).reshape(nnode_set,ndim**2,increments)
                e_hyd = np.einsum('ji',e_hyd).reshape(nnode_set,increments)
                e_VM = np.einsum('ji',e_VM).reshape(nnode_set,increments)
                s_VM = np.einsum('ji',s_VM).reshape(nnode_set,increments)
                FibreStress = np.einsum('lij',FibreStress).reshape(nnode_set,5,increments)
                collagen_density = collagen_density.reshape(nnode_set,increments)
            # Reshape tensors
            F = np.einsum('lijk',F).reshape(nnode_set,ndim**2,increments)
            J = J.reshape(nnode_set,increments)
            C = np.einsum('lijk',C).reshape(nnode_set,ndim**2,increments)
            detC = detC.reshape(nnode_set,increments)
            Cauchy = np.einsum('lijk',Cauchy).reshape(nnode_set,ndim**2,increments)
            p_hyd = np.einsum('ji',p_hyd).reshape(nnode_set,increments)

            if ndim == 2:
                C = C[:,[0,1,3],:]
                if self.gr_variables is None:
                    G = G[:,[0,1,3],:]
                else:
                    e = e[:,[0,1,3],:]
                Cauchy = Cauchy[:,[0,1,3],:]
            elif ndim == 3:
                C = C[:,[0,1,2,4,5,8],:]
                if self.gr_variables is None:
                    G = G[:,[0,1,2,4,5,8],:]
                else:
                    e = e[:,[0,1,2,4,5,8],:]
                Cauchy = Cauchy[:,[0,1,2,4,5,8],:]

            if fields == "mechanics" and ndim == 2 and self.gr_variables is None:
                augmented_sol[nstart:nend,:2,:]     = self.sol[self.materials[imat].node_set,:2,steps].reshape(augmented_sol[nstart:nend,:2,:].shape)
                augmented_sol[nstart:nend,2:6,:]    = F
                augmented_sol[nstart:nend,6:10,:]   = H
                augmented_sol[nstart:nend,10,:]     = J
                augmented_sol[nstart:nend,11:14,:]  = C
                augmented_sol[nstart:nend,14:17,:]  = G
                augmented_sol[nstart:nend,17,:]     = detC
                augmented_sol[nstart:nend,18:21,:]  = Cauchy
                augmented_sol[nstart:nend,21,:]     = p_hyd

            elif fields == "mechanics" and ndim == 3 and self.gr_variables is None:
                augmented_sol[nstart:nend,:3,:]     = self.sol[self.materials[imat].node_set,:3,steps].reshape(augmented_sol[nstart:nend,:3,:].shape)
                augmented_sol[nstart:nend,3:12,:]   = F
                augmented_sol[nstart:nend,12:21,:]  = H
                augmented_sol[nstart:nend,21,:]     = J
                augmented_sol[nstart:nend,22:28,:]  = C
                augmented_sol[nstart:nend,28:34,:]  = G
                augmented_sol[nstart:nend,34,:]     = detC
                augmented_sol[nstart:nend,35:41,:]  = Cauchy
                augmented_sol[nstart:nend,41,:]     = p_hyd
                augmented_sol[nstart:nend,42,:]     = s_VM

            elif fields == "mechanics" and ndim == 2 and self.gr_variables is not None:
                augmented_sol[nstart:nend,:2,:]     = self.sol[self.materials[imat].node_set,:2,steps].reshape(augmented_sol[nstart:nend,:2,:].shape)
                augmented_sol[nstart:nend,2:6,:]    = F
                augmented_sol[nstart:nend,6,:]      = J
                augmented_sol[nstart:nend,7:10,:]   = C
                augmented_sol[nstart:nend,10,:]     = detC
                augmented_sol[nstart:nend,11:14,:]  = e
                augmented_sol[nstart:nend,14,:]     = e_hyd
                augmented_sol[nstart:nend,15,:]     = e_VM
                augmented_sol[nstart:nend,16:19,:]  = Cauchy
                augmented_sol[nstart:nend,19,:]     = p_hyd
                augmented_sol[nstart:nend,20,:]     = s_VM
                if self.materials[imat].has_state_variables:
                    augmented_sol[nstart:nend,21:26,:]  = FibreStress
                    # Elastin and Muscle Densities
                    augmented_sol[nstart:nend,26:28,:]  = self.gr_variables[imat][:,5:7,steps].reshape(augmented_sol[nstart:nend,26:28,:].shape)
                    # Collagen Density
                    augmented_sol[nstart:nend,28,:]  = collagen_density
                    augmented_sol[nstart:nend,29,:]     = self.gr_variables[imat][:,11,steps].reshape(augmented_sol[nstart:nend,29,:].shape)   #Growth
                    augmented_sol[nstart:nend,30:35,:]  = self.gr_variables[imat][:,:5,steps].reshape(augmented_sol[nstart:nend,30:35,:].shape)   #Remodeling
                else:
                    augmented_sol[nstart:nend,21:26,:]  = np.zeros(FibreStress.shape,dtype=np.float64)
                    augmented_sol[nstart:nend,26:28,:]  = np.zeros(augmented_sol[nstart:nend,26:28,:].shape,dtype=np.float64) #EyMDensities
                    augmented_sol[nstart:nend,28,:]     = np.zeros(augmented_sol[nstart:nend,28,:].shape,dtype=np.float64)   #C Density
                    augmented_sol[nstart:nend,29,:]     = np.zeros(augmented_sol[nstart:nend,29,:].shape,dtype=np.float64)   #Growth
                    augmented_sol[nstart:nend,30:35,:]  = np.zeros(augmented_sol[nstart:nend,30:35,:].shape,dtype=np.float64)   #Remodeling

            elif fields == "mechanics" and ndim == 3 and self.gr_variables is not None:
                augmented_sol[nstart:nend,:3,:]     = self.sol[self.materials[imat].node_set,:3,steps].reshape(augmented_sol[nstart:nend,:3,:].shape)
                augmented_sol[nstart:nend,3:12,:]   = F
                augmented_sol[nstart:nend,12,:]     = J
                augmented_sol[nstart:nend,13:19,:]  = C
                augmented_sol[nstart:nend,19,:]     = detC
                augmented_sol[nstart:nend,20:26,:]  = e
                augmented_sol[nstart:nend,26,:]     = e_hyd
                augmented_sol[nstart:nend,27,:]     = e_VM
                augmented_sol[nstart:nend,28:34,:]  = Cauchy
                augmented_sol[nstart:nend,34,:]     = p_hyd
                augmented_sol[nstart:nend,35,:]     = s_VM
                if self.materials[imat].has_state_variables:
                    augmented_sol[nstart:nend,36:41,:]  = FibreStress
                    augmented_sol[nstart:nend,41:43,:]  = self.gr_variables[imat][:,5:7,steps].reshape(augmented_sol[nstart:nend,41:43,:].shape)
                    augmented_sol[nstart:nend,43,:]     = collagen_density
                    augmented_sol[nstart:nend,44,:]     = self.gr_variables[imat][:,11,steps].reshape(augmented_sol[nstart:nend,44,:].shape)
                    augmented_sol[nstart:nend,45:50,:]  = self.gr_variables[imat][:,:5,steps].reshape(augmented_sol[nstart:nend,45:50,:].shape)
                else:
                    augmented_sol[nstart:nend,36:41,:]  = np.zeros(FibreStress.shape,dtype=np.float64)
                    augmented_sol[nstart:nend,41:43,:]  = np.zeros(augmented_sol[nstart:nend,41:43,:].shape,dtype=np.float64)
                    augmented_sol[nstart:nend,43,:]     = np.zeros(augmented_sol[nstart:nend,43,:].shape,dtype=np.float64)
                    augmented_sol[nstart:nend,44,:]     = np.zeros(augmented_sol[nstart:nend,44,:].shape,dtype=np.float64)
                    augmented_sol[nstart:nend,45:50,:]  = np.zeros(augmented_sol[nstart:nend,45:50,:].shape,dtype=np.float64)

            nstart += nnode_set

        self.sol = augmented_sol

        return augmented_sol


    def WriteVTK(self,filename=None, quantity="all", configuration="deformed", steps=None, write_curved_mesh=True,
        interpolation_degree=10, ProjectionFlags=None, fmt="binary", equally_spaced=False, parallelise=False, time_problem=True):
        """Writes results to a VTK file for Paraview
            quantity = "all" means write all solution fields, otherwise specific quantities
            would be written based on augmented solution numbering order
            step - [list or np.1darray of sequentially aranged steps] which time steps/increments should be written
            inputs:
                fmt:                    [str] VTK writer format either "binary" or "xml".
                                        "xml" files do not support big vtk/vtu files
                                        typically greater than 2GB whereas "binary" files can.  Also "xml" writer is
                                        in-built whereas "binary" writer depends on evtk/pyevtk module
        """
        if fmt == "xml":
            pass
        elif fmt == "binary":
            try:
                from pyevtk.hl import pointsToVTK, linesToVTK, gridToVTK, unstructuredGridToVTK
                from pyevtk.vtk import VtkVertex, VtkLine, VtkTriangle, VtkQuad, VtkTetra, VtkPyramid, VtkHexahedron
            except ImportError:
                raise ImportError("Could not import evtk. Install it using 'pip install pyevtk'")
        else:
            raise ValueError("Writer format not understood")
        formatter = fmt

        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.sol is None:
            raise ValueError("solution not set for post-processing")

        if isinstance(quantity,int):
            if quantity>=self.sol.shape[1]:
                self.GetAugmentedSolution(parallelise=parallelise, time_problem=time_problem)
                if quantity >= self.sol.shape[1]:
                    raise ValueError('Plotting quantity not understood')
            iterator = range(quantity,quantity+1)
        elif isinstance(quantity,str):
            if quantity=="all":
                self.GetAugmentedSolution(parallelise=parallelise, time_problem=time_problem)
                iterator = range(self.sol.shape[1])
            else:
                raise ValueError('Plotting quantity not understood')
        elif isinstance(quantity,list):
            requires_augmented_solution = False
            for i in quantity:
                if i >= self.sol.shape[1]:
                    requires_augmented_solution = True
                    break
            if requires_augmented_solution:
                self.GetAugmentedSolution(parallelise=parallelise, time_problem=time_problem)
            iterator = quantity
        else:
            raise ValueError('Writing quantity not understood')

        if filename is None:
            warn("file name not specified. I am going to write in the current directory")
            filename = PWD(__file__) + "/output.vtu"
        elif filename is not None:
            if isinstance(filename,str) is False:
                raise ValueError("file name should be a string")
        if ".vtu" in filename and fmt is "binary":
            filename  = filename.split('.')[0]

        C = self.mesh.InferPolynomialDegree() - 1
        if C == 0:
            write_curved_mesh = False

        if self.sol.shape[1]>3:
            nnode = 0
            for imat in range(len(self.materials)):
                nnode += self.materials[imat].node_set.shape[0]
                self.materials[imat].ConnectivityOfMaterial(self.mesh)

            points = np.zeros((nnode,self.mesh.points.shape[1]),dtype=np.float64)
            elements = np.zeros((self.mesh.elements.shape[0],self.mesh.elements.shape[1]),dtype=np.int64)
            nstart = 0; nend = 0
            elstart = 0; elend = 0
            for imat in range(len(self.materials)):
                nend += self.materials[imat].node_set.shape[0]
                elend += self.materials[imat].element_set.shape[0]
                points[nstart:nend,:] = self.mesh.points[self.materials[imat].node_set]
                elements[elstart:elend,:] = nstart + self.materials[imat].elements
                nstart += self.materials[imat].node_set.shape[0]
                elstart += self.materials[imat].element_set.shape[0]
            self.mesh.points = points
            self.mesh.elements = elements

        # GET LINEAR MESH & SOLUTION
        lmesh, sol = self.mesh.GetLinearMesh(remap=True,solution=self.sol)

        if lmesh.element_type =='tri':
            cellflag = 5
            offset = 3
            actual_ndim = 2
        elif lmesh.element_type =='quad':
            cellflag = 9
            offset = 4
            actual_ndim = 2
        if lmesh.element_type =='tet':
            cellflag = 10
            offset = 4
            actual_ndim = 3
        elif lmesh.element_type == 'hex':
            cellflag = 12
            offset = 8
            actual_ndim = 3
        actual_ndim = lmesh.points.shape[1]

        ndim = lmesh.points.shape[1]
        q_names = [self.QuantityNamer(quant, print_name=True) for quant in iterator]

        TimeIncrement = self.sol.shape[2]

        increments = range(TimeIncrement)
        if steps is not None:
            increments = steps

        if write_curved_mesh is False:
            parallelise = False
        if len(increments) == 1:
            parallelise = False

        # PARALLEL MODE
        if parallelise:
            from multiprocessing import Pool, cpu_count
            from contextlib import closing
            if self.ncpu is None:
                self.ncpu = cpu_count()
            if self.parallel_model is None:
                self.parallel_model = "pool"

            increments = np.array(increments).flatten()
            partitioned_steps = np.array_split(increments,self.ncpu)

            pps = []
            for ip in range(len(partitioned_steps)):
                pp = deepcopy(self)
                # pp.sol = self.sol[:,:,partitioned_steps[ip]]
                pp.number_of_load_increments = len(partitioned_steps[ip])
                zipper_object = ParallelVTKWriterZipper(pp,
                    filename=filename, quantity=quantity, configuration=configuration,
                    steps=partitioned_steps[ip], write_curved_mesh=write_curved_mesh,
                    interpolation_degree=interpolation_degree,
                    ProjectionFlags=ProjectionFlags, fmt=fmt,
                    equally_spaced=equally_spaced)
                pps.append(zipper_object)

            # Pool
            if self.parallel_model == "pool":
                with closing(Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelWriteVTK,pps)
                    pool.terminate()

            # Thread Pool
            elif self.parallel_model == "thread_pool":
                import multiprocessing.dummy
                with closing(multiprocessing.dummy.Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelWriteVTK,pps)
                    pool.terminate()
            else:
                raise ValueError("Parallel model not understood")


            # # Serial
            # for ip in range(len(partitioned_steps)):
            #     ParallelWriteVTK(pps[ip])

            return

        if write_curved_mesh:

            if lmesh.element_type =='tet':
                cellflag = 5
                tmesh = PostProcess.TessellateTets(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_on_faces=False, plot_points=True,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags,
                    EquallySpacedPoints=equally_spaced)
            elif lmesh.element_type =='hex':
                cellflag = 5
                tmesh = PostProcess.TessellateHexes(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_on_faces=False, plot_points=True,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags,
                    EquallySpacedPoints=equally_spaced)
            elif lmesh.element_type =='quad':
                cellflag = 5
                tmesh = PostProcess.TessellateQuads(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_points=True, EquallySpacedPoints=equally_spaced,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags)
            elif lmesh.element_type =='tri':
                cellflag = 5
                tmesh = PostProcess.TessellateTris(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_points=True, EquallySpacedPoints=equally_spaced,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags)
            else:
                raise ValueError('Element type not understood')

            nsize = tmesh.nsize
            if hasattr(tmesh,'nface'):
                # FOR 3D ELEMENTS E.G. TETS AND HEXES
                nface = tmesh.nface
            else:
                tmesh.smesh = self.mesh
                tmesh.faces_to_plot = tmesh.smesh.elements
                nface = tmesh.smesh.elements.shape[0]
                tmesh.smesh.GetEdges()

                connections_elements = np.arange(tmesh.x_edges.size).reshape(tmesh.x_edges.shape[1],tmesh.x_edges.shape[0])
                connections = np.zeros((tmesh.x_edges.size,2),dtype=np.int64)
                for i in range(connections_elements.shape[0]):
                    connections[i*(tmesh.x_edges.shape[0]-1):(i+1)*(tmesh.x_edges.shape[0]-1),0] = connections_elements[i,:-1]
                    connections[i*(tmesh.x_edges.shape[0]-1):(i+1)*(tmesh.x_edges.shape[0]-1),1] = connections_elements[i,1:]
                connections = connections[:(i+1)*(tmesh.x_edges.shape[0]-1),:]
                tmesh.connections = connections

            un_faces_to_plot = np.unique(tmesh.faces_to_plot)
            fail_flag = False
            try:
                ssol = self.sol[un_faces_to_plot,:,:]
            except:
                fail_flag = True

            if fail_flag is False:
                if tmesh.smesh.elements.max() > un_faces_to_plot.shape[0]:
                    ssol = self.sol
                    fail_flag = True
                    warn("Something went wrong with mesh tessellation for VTK writer. I will proceed anyway")

            if tmesh.smesh.all_edges.shape[0] > tmesh.edge_elements.shape[0]:
                tmesh.smesh.all_edges = tmesh.edge_elements
                fail_flag = True
                warn("Something went wrong with mesh tessellation for VTK writer. I will proceed anyway")


            for Increment in increments:

                extrapolated_sol = np.zeros((tmesh.points.shape[0], self.sol.shape[1]))
                for ielem in range(nface):
                    extrapolated_sol[ielem*nsize:(ielem+1)*nsize,:] = np.dot(tmesh.bases_2,
                        ssol[tmesh.smesh.elements[ielem,:],:, Increment])

                if not fail_flag:
                    svpoints = self.mesh.points[np.unique(tmesh.faces_to_plot),:] + ssol[:,:tmesh.points.shape[1],Increment]
                else:
                    svpoints = self.mesh.points + ssol[:,:tmesh.points.shape[1],Increment]

                for iedge in range(tmesh.smesh.all_edges.shape[0]):
                    ielem = tmesh.edge_elements[iedge,0]
                    edge = tmesh.smesh.elements[ielem,tmesh.reference_edges[tmesh.edge_elements[iedge,1],:]]
                    coord_edge = svpoints[edge,:]
                    if tmesh.points.shape[1] == 3:
                        tmesh.x_edges[:,iedge], tmesh.y_edges[:,iedge], tmesh.z_edges[:,iedge] = np.dot(coord_edge.T,tmesh.bases_1)
                    elif tmesh.points.shape[1] == 2:
                        tmesh.x_edges[:,iedge], tmesh.y_edges[:,iedge] = np.dot(coord_edge.T,tmesh.bases_1)

                if tmesh.points.shape[1] == 3:
                    edge_coords = np.concatenate((tmesh.x_edges.T.copy().flatten()[:,None],
                        tmesh.y_edges.T.copy().flatten()[:,None],
                        tmesh.z_edges.T.copy().flatten()[:,None]),axis=1)
                elif tmesh.points.shape[1] == 2:
                    edge_coords = np.concatenate((tmesh.x_edges.T.copy().flatten()[:,None],
                        tmesh.y_edges.T.copy().flatten()[:,None], np.zeros_like(tmesh.y_edges.T.copy().flatten()[:,None])),axis=1)
                    svpoints = np.concatenate((svpoints, np.zeros((svpoints.shape[0],1))),axis=1)

                if formatter == "xml":

                    vtk_writer.write_vtu(Verts=edge_coords,
                        Cells={3:tmesh.connections},
                        fname=filename.split('.')[0]+'_curved_lines_increment_'+str(Increment)+'.vtu')

                    vtk_writer.write_vtu(Verts=svpoints,
                        Cells={1:np.arange(svpoints.shape[0])},
                        fname=filename.split('.')[0]+'_curved_points_increment_'+str(Increment)+'.vtu')

                    for quant in iterator:
                        vtk_writer.write_vtu(Verts=tmesh.points+extrapolated_sol[:,:ndim],
                            Cells={cellflag:tmesh.elements}, pdata=extrapolated_sol[:,quant],
                            fname=filename.split('.')[0]+'_curved_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')

                elif formatter == "binary":

                    unstructuredGridToVTK(filename.split('.')[0]+'_curved_lines_increment_'+str(Increment),
                        np.ascontiguousarray(edge_coords[:,0]), np.ascontiguousarray(edge_coords[:,1]), np.ascontiguousarray(edge_coords[:,2]),
                        np.ascontiguousarray(tmesh.connections.ravel()), np.arange(0,2*tmesh.connections.shape[0],2)+2,
                        np.ones(tmesh.connections.shape[0])*3)

                    pointsToVTK(filename.split('.')[0]+'_curved_points_increment_'+str(Increment),
                        np.ascontiguousarray(svpoints[:,0]), np.ascontiguousarray(svpoints[:,1]), np.ascontiguousarray(svpoints[:,2]),
                        data=None)

                    if tmesh.points.shape[1] == 2:
                        points = np.zeros((tmesh.points.shape[0],3))
                        points[:,:2] = tmesh.points+extrapolated_sol[:,:ndim]
                    else:
                        points = tmesh.points+extrapolated_sol[:,:ndim]
                    for counter, quant in enumerate(iterator):
                        unstructuredGridToVTK(filename.split('.')[0]+'_curved_quantity_'+str(quant)+'_increment_'+str(Increment),
                            np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                            np.ascontiguousarray(tmesh.elements.ravel()), np.arange(0,3*tmesh.elements.shape[0],3)+3,
                            np.ones(tmesh.elements.shape[0])*cellflag,
                            pointData={q_names[counter]: np.ascontiguousarray(extrapolated_sol[:,quant])})

        else:

            if configuration == "original":
                for Increment in increments:
                    if formatter == "xml":
                        for quant in iterator:
                            vtk_writer.write_vtu(Verts=lmesh.points,
                                Cells={cellflag:lmesh.elements}, pdata=sol[:,quant,Increment],
                                fname=filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')
                    elif formatter == "binary":
                        # points = lmesh.points
                        if lmesh.InferSpatialDimension() == 2:
                            points = np.zeros((lmesh.points.shape[0],3))
                            points[:,:2] = lmesh.points
                        else:
                            points = lmesh.points
                        for counter, quant in enumerate(iterator):
                            unstructuredGridToVTK(filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment),
                                np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                                np.ascontiguousarray(lmesh.elements.ravel()), np.arange(0,offset*lmesh.nelem,offset)+offset,
                                np.ones(lmesh.nelem)*cellflag,
                                pointData={q_names[counter]: np.ascontiguousarray(sol[:,quant,Increment])})

            elif configuration == "deformed":
                for Increment in increments:
                    if formatter == "xml":
                        for quant in iterator:
                            vtk_writer.write_vtu(Verts=lmesh.points+sol[:,:ndim,Increment],
                                Cells={cellflag:lmesh.elements}, pdata=sol[:,quant,Increment],
                                fname=filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')
                    elif formatter == "binary":
                        if lmesh.InferSpatialDimension() == 2:
                            points = np.zeros((lmesh.points.shape[0],3))
                            points[:,:2] = lmesh.points + sol[:,:ndim,Increment]
                        else:
                            points = lmesh.points + sol[:,:ndim,Increment]

                        for counter, quant in enumerate(iterator):
                            unstructuredGridToVTK(filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment),
                                np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                                np.ascontiguousarray(lmesh.elements.ravel()), np.arange(0,offset*lmesh.nelem,offset)+offset,
                                np.ones(lmesh.nelem)*cellflag,
                                pointData={q_names[counter]: np.ascontiguousarray(sol[:,quant,Increment])})

        return

    @staticmethod
    def TessellateHexes(mesh, TotalDisp, QuantityToPlot=None, plot_on_faces=True,
        ProjectionFlags=None, interpolation_degree=20, EquallySpacedPoints=False,
        plot_points=False, plot_edges=True, plot_surfaces=True):

        """High order curved hexahedral surfaces mesh plots, based on high order nodal FEM.
            The equally spaced FEM points do not work as good as the Fekete points
        """



        from Kuru.QuadratureRules import GaussLobattoPointsQuad
        from Kuru.QuadratureRules.NumericIntegrator import GaussLobattoQuadrature
        from Kuru.QuadratureRules.EquallySpacedPoints import EquallySpacedPoints as EquallySpacedPointsnD
        from Kuru.MeshGeneration.NodeArrangement import NodeArrangementQuad
        from Kuru.FunctionSpace import Quad
        from Kuru.FunctionSpace.OneDimensional.Line import LagrangeGaussLobatto, Lagrange

        from copy import deepcopy
        from scipy.spatial import Delaunay

        assert mesh.element_type == "hex"

        # SINCE THIS IS A 3D PLOT
        ndim=3

        C = interpolation_degree
        p = C+1
        nsize = int((C+2)**2)
        CActual = mesh.InferPolynomialDegree() - 1
        nsize_2 = int((CActual+2)**2)

        GaussLobattoPoints = GaussLobattoPointsQuad(C)
        if EquallySpacedPoints:
            GaussLobattoPoints = EquallySpacedPointsnD(ndim,C)

        # BUILD DELAUNAY TRIANGULATION OF REFERENCE ELEMENTS
        TrianglesFunc = Delaunay(GaussLobattoPoints)
        Triangles = TrianglesFunc.simplices.copy()

        # GET EQUALLY-SPACED/GAUSS-LOBATTO POINTS FOR THE EDGES
        GaussLobattoPointsOneD = GaussLobattoQuadrature(C+2)[0].flatten()
        if EquallySpacedPoints:
            GaussLobattoPointsOneD = EquallySpacedPointsnD(ndim-1, C).flatten()

        BasesQuad = np.zeros((nsize_2,GaussLobattoPoints.shape[0]),dtype=np.float64)
        hpBases = Quad.LagrangeGaussLobatto
        for i in range(GaussLobattoPoints.shape[0]):
            BasesQuad[:,i] = hpBases(CActual,GaussLobattoPoints[i,0],GaussLobattoPoints[i,1])[:,0]

        BasesOneD = np.zeros((CActual+2,GaussLobattoPointsOneD.shape[0]),dtype=np.float64)
        for i in range(GaussLobattoPointsOneD.shape[0]):
            BasesOneD[:,i] = LagrangeGaussLobatto(CActual,GaussLobattoPointsOneD[i])[0]

        # GET ONLY THE FACES WHICH NEED TO BE PLOTTED
        if ProjectionFlags is None:
            faces_to_plot_flag = np.ones(mesh.faces.shape[0])
        else:
            faces_to_plot_flag = ProjectionFlags.flatten()

        # CHECK IF ALL FACES NEED TO BE PLOTTED OR ONLY BOUNDARY FACES
        if faces_to_plot_flag.shape[0] > mesh.faces.shape[0]:
            # ALL FACES
            corr_faces = mesh.all_faces
            # FOR MAPPING DATA E.G. SCALED JACOBIAN FROM ELEMENTS TO FACES
            face_elements = mesh.GetElementsFaceNumberingHex()

        elif faces_to_plot_flag.shape[0] == mesh.faces.shape[0]:
            # ONLY BOUNDARY FACES
            corr_faces = mesh.faces
            # FOR MAPPING DATA E.G. SCALED JACOBIAN FROM ELEMENTS TO FACES
            face_elements = mesh.GetElementsWithBoundaryFacesHex()
        else:
            # raise ValueError("I do not understand what you want to plot")
            corr_faces = mesh.all_faces
            face_elements = mesh.GetElementsFaceNumberingHex()

        faces_to_plot = corr_faces[faces_to_plot_flag.flatten()==1,:]

        if QuantityToPlot is not None and plot_on_faces:
            quantity_to_plot = QuantityToPlot[face_elements[faces_to_plot_flag.flatten()==1,0]]

        # BUILD MESH OF SURFACE
        smesh = Mesh()
        smesh.element_type = "quad"
        # smesh.elements = np.copy(corr_faces)
        smesh.elements = np.copy(faces_to_plot)
        smesh.nelem = smesh.elements.shape[0]
        smesh.points = mesh.points[np.unique(smesh.elements),:]


        # MAP TO ORIGIN
        unique_elements, inv = np.unique(smesh.elements,return_inverse=True)
        mapper = np.arange(unique_elements.shape[0])
        smesh.elements = mapper[inv].reshape(smesh.elements.shape)

        smesh.GetBoundaryEdgesQuad()
        smesh.GetEdgesQuad()
        edge_elements = smesh.GetElementsEdgeNumberingQuad()



        # GET EDGE ORDERING IN THE REFERENCE ELEMENT
        reference_edges = NodeArrangementQuad(CActual)[0]
        reference_edges = np.concatenate((reference_edges,reference_edges[:,1,None]),axis=1)
        reference_edges = np.delete(reference_edges,1,1)

        # GET EULERIAN GEOMETRY
        if TotalDisp.ndim == 3:
            vpoints = mesh.points + TotalDisp[:,:ndim,-1]
        elif TotalDisp.ndim == 2:
            vpoints = mesh.points + TotalDisp[:,:ndim]
        else:
            raise AssertionError("mesh points and displacment arrays are incompatible")

        # svpoints = vpoints[np.unique(mesh.faces),:]
        svpoints = vpoints[np.unique(faces_to_plot),:]
        del vpoints
        gc.collect()

        if plot_edges:
            # GET X, Y & Z OF CURVED EDGES
            x_edges = np.zeros((C+2,smesh.all_edges.shape[0]))
            y_edges = np.zeros((C+2,smesh.all_edges.shape[0]))
            z_edges = np.zeros((C+2,smesh.all_edges.shape[0]))

            for iedge in range(smesh.all_edges.shape[0]):
                ielem = edge_elements[iedge,0]
                edge = smesh.elements[ielem,reference_edges[edge_elements[iedge,1],:]]
                coord_edge = svpoints[edge,:]
                x_edges[:,iedge], y_edges[:,iedge], z_edges[:,iedge] = np.dot(coord_edge.T,BasesOneD)


            # PLOT CURVED EDGES
            connections_elements = np.arange(x_edges.size).reshape(x_edges.shape[1],x_edges.shape[0])
            connections = np.zeros((x_edges.size,2),dtype=np.int64)
            for i in range(connections_elements.shape[0]):
                connections[i*(x_edges.shape[0]-1):(i+1)*(x_edges.shape[0]-1),0] = connections_elements[i,:-1]
                connections[i*(x_edges.shape[0]-1):(i+1)*(x_edges.shape[0]-1),1] = connections_elements[i,1:]
            connections = connections[:(i+1)*(x_edges.shape[0]-1),:]

        # CURVED SURFACES
        if plot_surfaces:

            nface = smesh.elements.shape[0]
            nnode = nsize*nface
            nelem = Triangles.shape[0]*nface

            Xplot = np.zeros((nnode,3),dtype=np.float64)
            Tplot = np.zeros((nelem,3),dtype=np.int64)

            # FOR CURVED ELEMENTS
            for ielem in range(nface):
                Xplot[ielem*nsize:(ielem+1)*nsize,:] = np.dot(BasesQuad.T, svpoints[smesh.elements[ielem,:],:])
                Tplot[ielem*TrianglesFunc.nsimplex:(ielem+1)*TrianglesFunc.nsimplex,:] = Triangles + ielem*nsize

            if QuantityToPlot is not None:
                Uplot = np.zeros(nnode,dtype=np.float64)
                if plot_on_faces:
                    # for ielem in range(nface):
                    #     Uplot[ielem*nsize:(ielem+1)*nsize] = quantity_to_plot[ielem]
                    # NEW APPROACH FOR CELL DATA - CHECK
                    Uplot = np.zeros(nelem,dtype=np.float64)
                    for ielem in range(nface):
                        Uplot[ielem*TrianglesFunc.nsimplex:(ielem+1)*TrianglesFunc.nsimplex] = quantity_to_plot[ielem]
                else:
                    # IF QUANTITY IS DEFINED ON NODES
                    quantity = QuantityToPlot[np.unique(faces_to_plot)]
                    for ielem in range(nface):
                        Uplot[ielem*nsize:(ielem+1)*nsize] = np.dot(BasesQuad.T, quantity[smesh.elements[ielem,:]])



        # THIS IS NOT A FLORENCE MESH COMPLIANT MESH
        tmesh = Mesh()
        tmesh.element_type = "tri"
        if plot_surfaces:
            tmesh.elements = Tplot
            tmesh.points = Xplot
            if QuantityToPlot is not None:
                tmesh.quantity = Uplot
            tmesh.nelem = nelem
            tmesh.nnode = nnode
            tmesh.nface = nface
        tmesh.nsize = nsize
        tmesh.bases_1 = BasesOneD
        tmesh.bases_2 = BasesQuad.T

        tmesh.smesh = smesh
        tmesh.faces_to_plot = faces_to_plot
        tmesh.svpoints = svpoints

        if plot_edges:
            tmesh.x_edges = x_edges
            tmesh.y_edges = y_edges
            tmesh.z_edges = z_edges
            tmesh.connections = connections
            tmesh.edge_elements = edge_elements
            tmesh.reference_edges = reference_edges

        return tmesh

