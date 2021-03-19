from __future__ import print_function
import gc, os, sys
import numpy as np
import scipy as sp
import numpy.linalg as la
import scipy.linalg as sla
from numpy.linalg import norm
from time import time
from copy import deepcopy
from warnings import warn
from time import time

from Kuru.FiniteElements.Assembly import Assemble #, AssembleExplicit
from Kuru.FiniteElements.LocalAssembly.KinematicMeasures import *
from Kuru.FiniteElements.LocalAssembly._KinematicMeasures_ import _KinematicMeasures_
from Kuru import Mesh
from Kuru import FunctionSpace, QuadratureRule

__all__ = ["GrowthRemodelingIntegrator"]


class GrowthRemodelingIntegrator(object):
    """Base class for structural time integerators
    """

    def __init__(self, gain, turnover, density_turnover="self", degradation_at_line=True,
        degradation_at_point=False, degradation_point=None, aging_only=False, monitoring_node=0, 
        damage_spread_space=0.010, damage_spread_time=40.0, damage_axis=0, **kwargs):

        self.HomeostaticStress = None
        self.gain = gain
        self.turnover = turnover
        self.density_turnover = density_turnover

        self.degradation_at_line = degradation_at_line
        self.degradation_at_point = degradation_at_point
        self.degradation_point = degradation_point
        self.aging_only = aging_only

        self.monitoring_node = monitoring_node

        if degradation_at_point or aging_only:
            self.degradation_at_line = False
        if degradation_point is None and degradation_at_point:
            self.degradation_point = [0.,0.,0.]

        self.damage_spread_space = damage_spread_space
        self.damage_spread_time = damage_spread_time
        self.damage_axis = damage_axis


    def HomeostaticDistortion(self, fem_solver, formulation, TotalDisp, Increment):
        """ check the distortion of homeostasis"""

        dmesh = Mesh()
        dmesh.points = TotalDisp[:,:formulation.ndim,Increment]
        dmesh_bounds = dmesh.Bounds
        distortion = 100.0*np.sqrt(dmesh_bounds[1,0]**2+dmesh_bounds[1,1]**2+\
                dmesh_bounds[1,2]**2)/0.010
        if distortion<5.0:
            print("The Distortion in Homeostasis is: {}".format(distortion))
        else:
            print("The Distortion in Homeostasis is: {}".format(distortion))
            sys.exit("Growth and Remodeling solver stop, distortion in Homeostasis is to big")


    def LogSave(self, fem_solver, formulation, TotalDisp, Increment, materials, FibreStress, gr_materials):

        if fem_solver.print_incremental_log:
            # find the set of the node under surveillance
            imat = -1
            for i in range(gr_materials.shape[0]):
                if self.monitoring_node in materials[gr_materials[i]].node_set:
                    imat = gr_materials[i]
                    imat0 = i
                    inode = np.where(materials[imat].node_set==self.monitoring_node)[0][0]
                    break
            if imat is -1:
                print("Set of the node is not recognized. I will use material 0 and its node 0.")
                imat = gr_materials[0]
                imat0 = 0
                inode = 0

            dmesh = Mesh()
            dmesh.points = TotalDisp[:,:formulation.ndim,Increment]
            dmesh_bounds = dmesh.Bounds
            print("\nMinimum and maximum incremental solution values at increment {} are \n".\
                    format(Increment),dmesh_bounds)

            print("\nGrowth and Remodeling properties at node, {}".format(self.monitoring_node))
            print("Densities: {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}".\
                    format(materials[imat].state_variables[inode,14],materials[imat].state_variables[inode,15],\
                    materials[imat].state_variables[inode,16],materials[imat].state_variables[inode,17],\
                    materials[imat].state_variables[inode,18],materials[imat].state_variables[inode,19]))
            print("Remodeling: {:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}".\
                    format(materials[imat].state_variables[inode,9],materials[imat].state_variables[inode,10],\
                    materials[imat].state_variables[inode,11],materials[imat].state_variables[inode,12],\
                    materials[imat].state_variables[inode,13]))
            print("Growth: {:6.3f}".format(materials[imat].state_variables[inode,20]))
            print("FibreStress: {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}".\
                    format(FibreStress[imat0][inode,0],FibreStress[imat0][inode,1],FibreStress[imat0][inode,2],\
                    FibreStress[imat0][inode,3],FibreStress[imat0][inode,4]))

        # SAVE INCREMENTAL SOLUTION IF ASKED FOR
        if fem_solver.save_incremental_solution:
            # FOR BIG MESHES
            if Increment % fem_solver.incremental_solution_save_frequency !=0:
                return
            from scipy.io import savemat
            filename = fem_solver.incremental_solution_filename
            if filename is not None:
                if ".mat" in filename:
                    filename = filename.split(".")[0]
                savemat(filename+"_"+str(Increment),
                        {'solution':TotalDisp[:,:,Increment]},do_compression=True)
            else:
                raise ValueError("No file name provided to save incremental solution")

    def NewtonRaphson(self, function_spaces, formulation, solver, Increment, K, NodalForces, 
            Residual, mesh, Eulerx, materials, boundary_condition, AppliedDirichletInc, fem_solver):

        Tolerance = fem_solver.newton_raphson_tolerance
        LoadIncrement = fem_solver.number_of_load_increments
        Iter = 0
        fem_solver.iterative_norm_history = []

        # APPLY INCREMENTAL DIRICHLET PER LOAD STEP (THIS IS INCREMENTAL NOT ACCUMULATIVE)
        IncDirichlet = boundary_condition.UpdateFixDoFs(AppliedDirichletInc,
            K.shape[0],formulation.nvar)
        # UPDATE EULERIAN COORDINATE
        Eulerx += IncDirichlet[:,:formulation.ndim]

        while fem_solver.norm_residual > Tolerance or Iter==0:
            # GET THE REDUCED SYSTEM OF EQUATIONS
            K_b, F_b = boundary_condition.GetReducedMatrices(K,Residual)[:2]

            # SOLVE THE SYSTEM
            sol = solver.Solve(K_b,-F_b)

            # GET ITERATIVE SOLUTION
            dU = boundary_condition.UpdateFreeDoFs(sol,K.shape[0],formulation.nvar)

            # UPDATE THE EULERIAN COMPONENTS
            # UPDATE THE GEOMETRY
            Eulerx += dU[:,:formulation.ndim]

            # RE-ASSEMBLE - COMPUTE STIFFNESS AND INTERNAL TRACTION FORCES
            K, TractionForces = Assemble(fem_solver, function_spaces, formulation, mesh, materials,
                boundary_condition, Eulerx)[:2]
            # COMPUTE ROBIN STIFFNESS AND FORCES (EXTERNAL)
            K, TractionForces = boundary_condition.ComputeRobinForces(mesh, materials, function_spaces,
                fem_solver, Eulerx, K, TractionForces)

            # FIND THE RESIDUAL
            Residual[boundary_condition.columns_in] = TractionForces[boundary_condition.columns_in] - \
                NodalForces[boundary_condition.columns_in]

            # SAVE THE NORM
            fem_solver.abs_norm_residual = la.norm(Residual[boundary_condition.columns_in])
            if Iter==0:
                fem_solver.NormForces = la.norm(Residual[boundary_condition.columns_in])
            fem_solver.norm_residual = np.abs(la.norm(Residual[boundary_condition.columns_in])/fem_solver.NormForces)

            # SAVE THE NORM
            fem_solver.NRConvergence['Increment_'+str(Increment)] = np.append(fem_solver.NRConvergence[\
                    'Increment_'+str(Increment)],fem_solver.norm_residual)

            print("Iteration {} for increment {}.".format(Iter, Increment) +\
                " Residual (abs) {0:>16.7g}".format(fem_solver.abs_norm_residual),
                "\t Residual (rel) {0:>16.7g}".format(fem_solver.norm_residual))

            # BREAK BASED ON RELATIVE NORM
            if np.abs(fem_solver.abs_norm_residual) < Tolerance:
                break

            # BREAK BASED ON INCREMENTAL SOLUTION - KEEP IT AFTER UPDATE
            if norm(dU) <=  fem_solver.newton_raphson_solution_tolerance and Iter!=0:
                print("Incremental solution within tolerance i.e. norm(dU): {}".format(norm(dU)))
                break

            # UPDATE ITERATION NUMBER
            Iter +=1

            if Iter==fem_solver.maximum_iteration_for_newton_raphson:
                fem_solver.newton_raphson_failed_to_converge = True
                break
            if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
                fem_solver.newton_raphson_failed_to_converge = True
                break

            # IF BREAK WHEN NEWTON RAPHSON STAGNATES IS ACTIVATED
            if fem_solver.break_at_stagnation:
                fem_solver.iterative_norm_history.append(fem_solver.norm_residual)
                if Iter >= 5:
                    if np.mean(fem_solver.iterative_norm_history) < 1.:
                        break

            # USER DEFINED CRITERIA TO BREAK OUT OF NEWTON-RAPHSON
            if fem_solver.user_defined_break_func != None:
                if fem_solver.user_defined_break_func(Increment,Iter,fem_solver.norm_residual,fem_solver.abs_norm_residual, Tolerance):
                    break

            # USER DEFINED CRITERIA TO STOP NEWTON-RAPHSON AND THE WHOLE ANALYSIS
            if fem_solver.user_defined_stop_func != None:
                if fem_solver.user_defined_stop_func(Increment,Iter,fem_solver.norm_residual,fem_solver.abs_norm_residual, Tolerance):
                    fem_solver.newton_raphson_failed_to_converge = True
                    break


        return Eulerx, K, Residual

    def GetFibreStressAndSoftness(self, mesh, formulation, material, fem_solver, Eulerx, average_derived_quantities=True):

        """
            steps:          [list,np.1darray] for which time steps/increments the data should
                            be recovered
        """

        det = np.linalg.det
        inv = np.linalg.inv

        # GET THE UNDERLYING LINEAR MESH
        # lmesh = mesh.GetLinearMesh()
        C = mesh.InferPolynomialDegree() - 1
        ndim = mesh.InferSpatialDimension()

        nelem = mesh.elements.shape[0]; npoint = mesh.points.shape[0]
        nodeperelem = mesh.elements.shape[1]

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

        F = np.zeros((material.element_set.shape[0],nodeperelem,ndim,ndim))
        # DEFINE CONSTITUENT STRESSES FOR GROWTH-REMODELING PROBLEM
        ElemFibreStress = np.zeros((material.element_set.shape[0],nodeperelem,5))  # 5-fibres
        ElemSoftness = np.zeros((material.element_set.shape[0],nodeperelem,5))  # 5-fibres

        FibreStress = np.zeros((material.node_set.shape[0],5))
        Softness = np.zeros((material.node_set.shape[0],5))

        # LOOP OVER ELEMENTS
        for ielem in range(material.element_set.shape[0]):
            elem = material.element_set[ielem]
            # GET THE FIELDS AT THE ELEMENT LEVEL
            LagrangeElemCoords = mesh.points[mesh.elements[elem,:],:]
            EulerElemCoords = Eulerx[mesh.elements[elem,:],:]
            # GROWTH-REMODELING VALUES FOR THIS ELEMENT
            material.MappingStateVariables(mesh,Domain,elem)

            if material.has_low_level_dispatcher:
                # GET LOCAL KINEMATICS
                SpatialGradient, F[ielem,:,:,:], detJ, dV = _KinematicMeasures_(Jm, AllGauss[:,0],
                        LagrangeElemCoords, EulerElemCoords, requires_geometry_update)
                # PARAMETERS FOR INCOMPRESSIBILITY (MEAN DILATATION METHOD HU-WASHIZU)
                if material.is_incompressible:
                    MaterialVolume = np.sum(dV)
                    if fem_solver.has_growth_remodeling:
                        dve = np.true_divide(detJ,material.StateVariables[:,20])
                        CurrentVolume = np.sum(dve)
                    else:
                        CurrentVolume = np.sum(detJ)
                    material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                # COMPUTE FIBRE STRESS AND SOFTNESS
                ElemFibreStress[ielem,:,:],ElemSoftness[ielem,:,:] = material._ConstituentMeasures_(F[ielem,:,:,:],elem)

            else:
                # GAUSS LOOP IN VECTORISED FORM
                ParentGradientX = np.einsum('ijk,jl->kil', Jm, LagrangeElemCoords)
                # MATERIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla_0 (N)]
                MaterialGradient = np.einsum('ijk,kli->ijl', inv(ParentGradientX), Jm)
                # DEFORMATION GRADIENT TENSOR [\vec{x} \otimes \nabla_0 (N)]
                F[ielem,:,:,:] = np.einsum('ij,kli->kjl', EulerElemCoords, MaterialGradient)
                # COMPUTE REMAINING KINEMATIC MEASURES
                StrainTensors = KinematicMeasures(F[ielem,:,:,:], fem_solver.analysis_nature)

                # GEOMETRY UPDATE IS A MUST
                # MAPPING TENSOR [\partial\vec{X}/ \partial\vec{\varepsilon} (ndim x ndim)]
                ParentGradientx = np.einsum('ijk,jl->kil',Jm,EulerElemCoords)
                # SPATIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla (N)]
                SpatialGradient = np.einsum('ijk,kli->ilj',inv(ParentGradientx),Jm)
                # COMPUTE ONCE detJ (GOOD SPEEDUP COMPARED TO COMPUTING TWICE)
                detJ = np.einsum('i,i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)),
                    np.abs(StrainTensors['J']))

                # COMPUTE PARAMETERS FOR MEAN DILATATION METHOD, IT NEEDS TO BE BEFORE COMPUTE HESSIAN AND STRESS
                if material.is_incompressible:
                    dV = np.einsum('i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)))
                    MaterialVolume = np.sum(dV)
                    if fem_solver.has_growth_remodeling:
                        dve = np.true_divide(detJ,material.StateVariables[:,20])
                        CurrentVolume = np.sum(dve)
                    else:
                        CurrentVolume = np.sum(detJ)
                    material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                # LOOP OVER GAUSS POINTS
                for counter in range(AllGauss.shape[0]):
                    # COMPUTE FIBRE STRESS AND SOFTNESS
                    ElemFibreStress[ielem,counter,:],ElemSoftness[ielem,counter,:] = material.ConstituentMeasures(
                        StrainTensors,elem,counter)


        # COMPUTE THE COMMON/NEIGHBOUR NODES ONCE
        Elss, Poss = material.GetNodeCommonality()[:2]
        for inode in range(material.node_set.shape[0]):
            Els, Pos = Elss[inode], Poss[inode]
            ncommon_nodes = Els.shape[0]
            for uelem in range(ncommon_nodes):
                FibreStress[inode,:] += ElemFibreStress[Els[uelem],Pos[uelem],:]
                Softness[inode,:] += ElemSoftness[Els[uelem],Pos[uelem],:]

            # AVERAGE OUT
            FibreStress[inode,:] /= ncommon_nodes
            Softness[inode,:] /= ncommon_nodes


        return FibreStress,Softness

    def RatesGrowthRemodeling(self, mesh, material, FibreStress, Softness, imat):
        """ This are the rates of Growth and Rmodeling
        """

        if self.HomeostaticStress is None:
            raise ValueError("Homeostatic Stress is not fixed")

        #k_s = np.zeros((),dytpe=)
        k_s = self.gain/self.turnover
        Rates = np.zeros((material.node_set.shape[0],10),dtype=np.float64)

        # choose a mode of collagen addition, either self-driven or muscle-driven
        # each fibre density is driven by its own stress
        if self.density_turnover is "self":
            for node in range(material.node_set.shape[0]):
                for fibre in range(5):
                    if self.HomeostaticStress[imat][node,fibre] == 0.:
                        continue
                    DeltaStress = FibreStress[imat][node,fibre] - self.HomeostaticStress[imat][node,fibre]
                    # Fibre density rate
                    Rates[node,fibre+5] = k_s*material.state_variables[node,fibre+15]*\
                        DeltaStress/self.HomeostaticStress[imat][node,fibre]
                    # Fibre remodeling rate
                    Rates[node,fibre] = (Rates[node,fibre+5]/material.state_variables[node,fibre+15] + \
                        1./self.turnover)*DeltaStress*Softness[imat][node,fibre]

        # each fibre density is driven by muscle stress
        elif self.density_turnover is "muscle":
            for node in range(material.node_set.shape[0]):
                for fibre in range(5):
                    if self.HomeostaticStress[imat][node,fibre] == 0.:
                        continue
                    # assuming the material 0 is the media
                    DeltaStress_m = FibreStress[0][node,0] - self.HomeostaticStress[0][node,0]
                    DeltaStress = FibreStress[imat][node,fibre] - self.HomeostaticStress[imat][node,fibre]
                    # Fibre density rate
                    Rates[node,fibre+5] = k_s*material.state_variables[node,fibre+15]*\
                        DeltaStress_m/self.HomeostaticStress[0][node,0]
                    # Fibre remodeling rate
                    Rates[node,fibre] = (Rates[node,fibre+5]/material.state_variables[node,fibre+15] + \
                        1./self.turnover)*DeltaStress*Softness[imat][node,fibre]

        # each fibre density is driven by its own stress and mass is added just when de delta is positive
        elif self.density_turnover is "self_sgn":
            for node in range(material.node_set.shape[0]):
                for fibre in range(5):
                    DeltaStress = FibreStress[imat][node,fibre] - self.HomeostaticStress[imat][node,fibre]
                    if DeltaStress > 0.0:
                        # Fibre density rate
                        Rates[node,fibre+5] = k_s*material.state_variables[node,fibre+15]*\
                            DeltaStress/self.HomeostaticStress[imat][node,fibre]
                        # Fibre remodeling rate
                        Rates[node,fibre] = (Rates[node,fibre+5]/material.state_variables[node,fibre+15] + \
                            1./self.turnover)*DeltaStress*Softness[imat][node,fibre]
                    else:
                        # Fibre density rate
                        Rates[node,fibre+5] = 0.0
                        # Fibre remodeling rate
                        Rates[node,fibre] = 0.0


        return Rates


