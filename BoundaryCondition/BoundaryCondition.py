from __future__ import print_function
import sys
import numpy as np #, scipy as sp, os, gc
from copy import deepcopy
#from warnings import warn
from time import time

class BoundaryCondition(object):
    """Base class for applying all types of boundary conditions"""

    def __init__(self,
        surface_identification_algorithm='minimisation',
        modify_linear_mesh_on_projection=False,
        project_on_curves=True,
        activate_bounding_box=False,
        bounding_box_padding=1e-3,
        has_planar_surfaces=True,
        solve_for_planar_faces=True,
        save_dirichlet_data=False,
        save_nurbs_data=False,
        filename=None,
        read_dirichlet_from_file=False,
        make_loading="ramp",
        compound_dirichlet_bcs=False
        ):

        # TYPE OF BOUNDARY: straight or nurbs
        self.boundary_type = 'straight'
        self.dirichlet_data_applied_at = 'node' # or 'faces'
        self.neumann_data_applied_at = 'node' # or 'faces'
        self.requires_cad = False
        self.cad_file = None
        # PROJECTION TYPE FOR CAD EITHER orthogonal OR arc_length
        self.projection_type = 'orthogonal'
        # WHAT TYPE OF ARC LENGTH BASED PROJECTION, EITHER 'equal' OR 'fekete'
        self.nodal_spacing_for_cad = 'equal'
        self.project_on_curves = project_on_curves
        self.scale_mesh_on_projection = False
        self.scale_value_on_projection = 1.0
        self.condition_for_projection = 1.0e20
        self.has_planar_surfaces = False
        self.solve_for_planar_faces = solve_for_planar_faces
        self.projection_flags = None
        # FIX DEGREES OF FREEDOM EVERY WHERE CAD PROJECTION IS NOT APPLIED
        self.fix_dof_elsewhere = True
        # FOR 3D ARC-LENGTH PROJECTION
        self.orthogonal_fallback_tolerance = 1.0
        # WHICH ALGORITHM TO USE FOR SURFACE IDENTIFICATION, EITHER 'minimisation' or 'pure_projection'
        self.surface_identification_algorithm = surface_identification_algorithm
        # MODIFY LINEAR MESH ON PROJECTION
        self.modify_linear_mesh_on_projection = modify_linear_mesh_on_projection
        # COMPUTE A BOUNDING BOX FOR EACH CAD SURFACE
        self.activate_bounding_box = activate_bounding_box
        self.bounding_box_padding = float(bounding_box_padding)

        # FOR IGAKit WRAPPER
        self.nurbs_info = None
        self.nurbs_condition = None

        self.analysis_type = 'static'
        self.analysis_nature = 'linear'

        self.dirichlet_flags = None
        self.applied_dirichlet = None
        self.is_dirichlet_computed = False
        self.columns_out = None
        self.columns_in = None
        self.save_dirichlet_data = save_dirichlet_data
        self.save_nurbs_data = save_nurbs_data
        self.filename = filename
        self.read_dirichlet_from_file = read_dirichlet_from_file

        self.neumann_flags = None
        self.applied_neumann = None
        self.is_applied_neumann_shape_functions_computed = False

        self.pressure_flags = None
        self.applied_pressure = None
        self.pressure_increment = 1.0
        self.spring_flags = None
        self.applied_spring = None
        self.connector_elements = None
        self.applied_connector = None
        self.is_body_force_shape_functions_computed = False

        self.make_loading = make_loading # "ramp" or "constant"
        self.has_step_wise_dirichlet_loading = False
        self.step_wise_dirichlet_data = None
        self.has_step_wise_neumann_loading = False
        self.step_wise_neumann_data = None

        self.compound_dirichlet_bcs = compound_dirichlet_bcs

        # STORE A COPY OF SELF AT THE START TO RESET TO AT THE END
        self.__save_state__()
        # FOR INTERNAL PURPOSES WHEN WE DO NOT WANT TO REST
        self.do_not_reset = True

    def __save_state__(self):
        self.__initialdict__ = deepcopy(self.__dict__)

    def SetDirichletCriteria(self, func, *args, **kwargs):
        """Applies user defined Dirichlet data to self
        """

        if "apply" in kwargs.keys():
            del kwargs["apply"]
            self.has_step_wise_dirichlet_loading = True
            self.step_wise_dirichlet_data = {'func':func, 'args': args, 'kwargs': kwargs}
            self.dirichlet_flags = func(0, *args, **kwargs)
            return self.dirichlet_flags

        self.dirichlet_flags = func(*args, **kwargs)
        return self.dirichlet_flags

    def SetNeumannCriteria(self, func, *args, **kwargs):
        """Applies user defined Neumann data to self
        """

        if "apply" in kwargs.keys():
            del kwargs["apply"]
            self.has_step_wise_neumann_loading = True
            self.step_wise_neumann_data = {'func':func, 'args': args, 'kwargs': kwargs}
            tups = func(0, *args, **kwargs)
        else:
            tups = func(*args, **kwargs)

        if not isinstance(tups,tuple) and self.neumann_data_applied_at == "node":
            self.neumann_flags = tups
            return self.neumann_flags
        else:
            self.neumann_data_applied_at == "face"
            if len(tups) !=2:
                raise ValueError("User-defined Neumann criterion function {} "
                    "should return one flag and one data array".format(func.__name__))
            self.neumann_flags = tups[0]
            self.applied_neumann = tups[1]
            return tups

    def SetRobinCriteria(self, func, *args, **kwargs):
        """Applies user defined Robin data to self, just working on surfaces
        """

        dics = func(*args, **kwargs)

        if isinstance(dics,dict):
            self.RobinLoadSelector(dics)
        elif isinstance(dics,tuple):
            for idic in range(len(dics)):
                if isinstance(dics[idic],dict):
                    self.RobinLoadSelector(dics[idic])
                else:
                    raise ValueError("User-defined Robin criterion function {} "
                        "should return dictionary or tuple(dict,dict,...)".format(func.__name__))
        else:
            raise ValueError("User-defined Robin criterion function {} "
                "should return dictionary or tuple".format(func.__name__))

        return dics

    def RobinLoadSelector(self, tups):
        if tups['type'] == 'Pressure':
            self.pressure_flags = tups['flags']
            self.applied_pressure = tups['data']
        elif tups['type'] == 'Spring':
            self.spring_flags = tups['flags']
            self.applied_spring = tups['data']
        elif tups['type'] == 'Connector':
            self.connector_elements = tups['elements']
            self.applied_connector = tups['data']
        elif tups['type'] == 'Dashpot':
            raise ValueError("Surrounding viscoelastic effects not implemented yet")
        else:
            raise ValueError("Type force {} not understood or not available. "
                "Types are Pressure, Spring, SpringJoint and Dashpot.".format(tups['type']))


    def GetDirichletBoundaryConditions(self, formulation, mesh, materials=None, solver=None, fem_solver=None):

        nvar = formulation.nvar
        ndim = formulation.ndim
        self.columns_in, self.applied_dirichlet = [], []

        #----------------------------------------------------------------------------------------------------#
        #-------------------------------------- NURBS BASED SOLUTION ----------------------------------------#
        #----------------------------------------------------------------------------------------------------#
        if self.boundary_type == 'nurbs':

            tCAD = time()

            if self.read_dirichlet_from_file is False:

                if not self.is_dirichlet_computed:
                    # GET DIRICHLET BOUNDARY CONDITIONS BASED ON THE EXACT GEOMETRY FROM CAD
                    if self.requires_cad:
                        # CALL POSTMESH WRAPPER
                        nodesDBC, Dirichlet = self.PostMeshWrapper(formulation, mesh, materials, solver, fem_solver)
                else:
                    nodesDBC, Dirichlet = self.nodesDBC, self.Dirichlet


                # GET DIRICHLET DoFs
                self.columns_out = (np.repeat(nodesDBC,nvar,axis=1)*nvar +\
                 np.tile(np.arange(nvar)[None,:],nodesDBC.shape[0]).reshape(nodesDBC.shape[0],formulation.ndim)).ravel()
                self.applied_dirichlet = Dirichlet.ravel()


                # FIX THE DOF IN THE REST OF THE BOUNDARY
                if self.fix_dof_elsewhere:
                    if ndim==2:
                        rest_dofs = np.setdiff1d(np.unique(mesh.edges),nodesDBC)
                    elif ndim==3:
                        rest_dofs = np.setdiff1d(np.unique(mesh.faces),nodesDBC)

                    rest_out = np.repeat(rest_dofs,nvar)*nvar + np.tile(np.arange(nvar),rest_dofs.shape[0])
                    rest_app = np.zeros(rest_dofs.shape[0]*nvar)

                    self.columns_out = np.concatenate((self.columns_out,rest_out)).astype(np.int64)
                    self.applied_dirichlet = np.concatenate((self.applied_dirichlet,rest_app))


                print('Finished identifying Dirichlet boundary conditions from CAD geometry.',
                    ' Time taken', time()-tCAD, 'seconds')

            else:

                end = -3
                self.applied_dirichlet = np.loadtxt(mesh.filename.split(".")[0][:end]+"_dirichlet.dat",  dtype=np.float64)
                self.columns_out = np.loadtxt(mesh.filename.split(".")[0][:end]+"_columns_out.dat")

                print('Finished identifying Dirichlet boundary conditions from CAD geometry.',
                    ' Time taken', time()-tCAD, 'seconds')

        #----------------------------------------------------------------------------------------------------#
        #------------------------------------- NON-NURBS BASED SOLUTION -------------------------------------#
        #----------------------------------------------------------------------------------------------------#

        elif self.boundary_type == 'straight' or self.boundary_type == 'mixed':
            # IF DIRICHLET BOUNDARY CONDITIONS ARE APPLIED DIRECTLY AT NODES
            if self.dirichlet_flags is None:
                raise RuntimeError("Dirichlet boundary conditions are not set for the analysis")

            if self.dirichlet_data_applied_at == 'node':
                if self.analysis_type == "dynamic":
                    # FOR DYNAMIC ANALYSIS IT IS ASSUMED THAT
                    # self.columns_in and self.columns_out DO NOT CHANGE
                    # DURING THE ANALYSIS
                    if self.dirichlet_flags.ndim == 3:
                        flat_dirich = self.dirichlet_flags[:,:,0].ravel()
                        self.columns_out = np.arange(self.dirichlet_flags[:,:,0].size)[~np.isnan(flat_dirich)]
                        self.applied_dirichlet = np.zeros((self.columns_out.shape[0],self.dirichlet_flags.shape[2]))

                        for step in range(self.dirichlet_flags.shape[2]):
                            flat_dirich = self.dirichlet_flags[:,:,step].ravel()
                            self.applied_dirichlet[:,step] = flat_dirich[~np.isnan(flat_dirich)]

                    elif self.dirichlet_flags.ndim == 2:
                        flat_dirich = self.dirichlet_flags.ravel()
                        self.columns_out = np.arange(self.dirichlet_flags.size)[~np.isnan(flat_dirich)]
                        self.applied_dirichlet = flat_dirich[~np.isnan(flat_dirich)]
                    else:
                        raise ValueError("Incorrect Dirichlet flags for dynamic analysis")

                else:
                    flat_dirich = self.dirichlet_flags.ravel()
                    self.columns_out = np.arange(self.dirichlet_flags.size)[~np.isnan(flat_dirich)]
                    self.applied_dirichlet = flat_dirich[~np.isnan(flat_dirich)]

        # GENERAL PROCEDURE - GET REDUCED MATRICES FOR FINAL SOLUTION
        self.columns_out = self.columns_out.astype(np.int64)
        self.columns_in = np.delete(np.arange(0,nvar*mesh.points.shape[0]),self.columns_out)

        if self.columns_in.shape[0] == 0:
            warn("No Dirichlet boundary conditions have been applied. The system is unconstrained")
        if self.columns_out.shape[0] == 0:
            warn("Dirichlet boundary conditions have been applied on the entire mesh")

        if self.save_dirichlet_data:
            from scipy.io import savemat
            diri_dict = {'columns_in':self.columns_in,
                'columns_out':self.columns_out,
                'applied_dirichlet':self.applied_dirichlet}
            savemat(self.filename,diri_dict, do_compression=True)


    def ComputeNeumannForces(self, mesh, materials, function_spaces, compute_traction_forces=True, compute_body_forces=False):
        """Compute/assemble traction and body forces"""

        if self.neumann_flags is None:
            return np.zeros((mesh.points.shape[0]*materials[0].nvar,1),dtype=np.float64)

        nvar = materials[0].nvar
        ndim = mesh.InferSpatialDimension()

        if self.neumann_flags.shape[0] == mesh.points.shape[0]:
            self.neumann_data_applied_at = "node"
        else:
            if ndim==3:
                if self.neumann_flags.shape[0] == mesh.faces.shape[0]:
                    self.neumann_data_applied_at = "face"
            elif ndim==2:
                if self.neumann_flags.shape[0] == mesh.edges.shape[0]:
                    self.neumann_data_applied_at = "face"


        if self.neumann_data_applied_at == 'face':
            from Kuru.FiniteElements.Assembly import AssembleForces
            if not isinstance(function_spaces,tuple):
                raise ValueError("Boundary functional spaces not available for computing Neumman and body forces")
            else:
                # CHECK IF A FUNCTION SPACE FOR BOUNDARY EXISTS - SAFEGAURDS AGAINST FORMULATIONS THAT DO NO PROVIDE ONE
                has_boundary_spaces = False
                for fs in function_spaces:
                    if ndim == 3 and fs.ndim == 2:
                        has_boundary_spaces = True
                        break
                    elif ndim == 2 and fs.ndim == 1:
                        has_boundary_spaces = True
                        break
                if not has_boundary_spaces:
                    from Kuru import QuadratureRule, FunctionSpace
                    # COMPUTE BOUNDARY FUNCTIONAL SPACES
                    p = mesh.InferPolynomialDegree()
                    bquadrature = QuadratureRule(optimal=3, norder=2*p+1,
                        mesh_type=mesh.boundary_element_type, is_flattened=False)
                    bfunction_space = FunctionSpace(mesh.CreateDummyLowerDimensionalMesh(),
                        bquadrature, p=p, equally_spaced=mesh.IsEquallySpaced, use_optimal_quadrature=False)
                    function_spaces = (function_spaces[0],bfunction_space)
                    # raise ValueError("Boundary functional spaces not available for computing Neumman and body forces")

            t_tassembly = time()
            if self.analysis_type == "static":
                F = AssembleForces(self, mesh, materials, function_spaces,
                    compute_traction_forces=compute_traction_forces, compute_body_forces=compute_body_forces)
            elif self.analysis_type == "dynamic":
                if self.neumann_flags.ndim==2:
                    # THE POSITION OF NEUMANN DATA APPLIED AT FACES CAN CHANGE DYNAMICALLY
                    tmp_flags = np.copy(self.neumann_flags)
                    tmp_data = np.copy(self.applied_neumann)
                    F = np.zeros((mesh.points.shape[0]*nvar,self.neumann_flags.shape[1]))
                    for step in range(self.neumann_flags.shape[1]):
                        self.neumann_flags = tmp_flags[:,step]
                        self.applied_neumann = tmp_data[:,:,step]
                        F[:,step] = AssembleForces(self, mesh, materials, function_spaces,
                            compute_traction_forces=compute_traction_forces, compute_body_forces=compute_body_forces).flatten()

                    self.neumann_flags = tmp_flags
                    self.applied_neumann = tmp_data
                else:
                    # THE POSITION OF NEUMANN DATA APPLIED AT FACES CAN CHANGE DYNAMICALLY
                    F = AssembleForces(self, mesh, materials, function_spaces,
                            compute_traction_forces=compute_traction_forces, compute_body_forces=compute_body_forces).flatten()

            print("Assembled external traction forces. Time elapsed is {} seconds".format(time()-t_tassembly))


        elif self.neumann_data_applied_at == 'node':
            # A DIRICHLET TYPE METHODOLGY FOR APPLYING NEUMANN BOUNDARY CONDITONS (i.e. AT NODES)
            if self.analysis_type == "dynamic":
                if self.neumann_flags.ndim ==3:
                    # FOR DYNAMIC ANALYSIS IT IS ASSUMED THAT
                    # to_apply DOOES NOT CHANGE DURING THE ANALYSIS
                    flat_neu = self.neumann_flags[:,:,0].ravel()
                    to_apply = np.arange(self.neumann_flags[:,:,0].size)[~np.isnan(flat_neu)]
                    F = np.zeros((mesh.points.shape[0]*nvar,self.neumann_flags.shape[2]))

                    for step in range(self.neumann_flags.shape[2]):
                        flat_neu = self.neumann_flags[:,:,step].ravel()
                        to_apply = np.arange(self.neumann_flags[:,:,step].size)[~np.isnan(flat_neu)]
                        F[to_apply,step] = flat_neu[~np.isnan(flat_neu)]
                else:
                    F = np.zeros((mesh.points.shape[0]*nvar,1))
                    flat_neu = self.neumann_flags.ravel()
                    to_apply = np.arange(self.neumann_flags.size)[~np.isnan(flat_neu)]
                    applied_neumann = flat_neu[~np.isnan(flat_neu)]
                    F[to_apply,0] = applied_neumann
            else:
                F = np.zeros((mesh.points.shape[0]*nvar,1))
                flat_neu = self.neumann_flags.ravel()
                to_apply = np.arange(self.neumann_flags.size)[~np.isnan(flat_neu)]
                applied_neumann = flat_neu[~np.isnan(flat_neu)]
                F[to_apply,0] = applied_neumann

        return F

    def ComputeRobinForces(self, mesh, materials, function_spaces, fem_solver, Eulerx, stiffness, F):
        """Compute/assemble traction and body forces"""

        from Kuru.FiniteElements.Assembly import AssembleRobinForces
        if not self.pressure_flags is None:
            K_pressure, F_pressure = AssembleRobinForces(self, mesh,
                materials[0], function_spaces, fem_solver, Eulerx, 'pressure')
            stiffness -= K_pressure
            F -= F_pressure[:,None]
        if not self.spring_flags is None:
            K_spring, F_spring = AssembleRobinForces(self, mesh,
                materials[0], function_spaces, fem_solver, Eulerx, 'spring')
            stiffness += K_spring
            F += F_spring[:,None]
        if not self.connector_elements is None:
            K_connector, F_connector = AssembleRobinForces(self, mesh,
                materials[0], function_spaces, fem_solver, Eulerx, 'connector')
            stiffness += K_connector
            F += F_connector[:,None]

        return stiffness, F

    def GetReducedMatrices(self, stiffness, F, mass=None, only_residual=False):

        # GET REDUCED FORCE VECTOR
        F_b = F[self.columns_in,0]
        if only_residual:
            return F_b

        # GET REDUCED STIFFNESS MATRIX
        stiffness_b = stiffness[self.columns_in,:][:,self.columns_in]

        # GET REDUCED MASS MATRIX
        mass_b = np.array([])

        return stiffness_b, F_b, mass_b

    def ApplyDirichletGetReducedMatrices(self, stiffness, F, AppliedDirichlet, LoadFactor=1., mass=None, only_residual=False):
        """AppliedDirichlet is a non-member because it can be external incremental Dirichlet,
            which is currently not implemented as member of BoundaryCondition. F also does not
            correspond to Dirichlet forces, as it can be residual in incrementally linearised
            framework.
        """

        # # APPLY DIRICHLET BOUNDARY CONDITIONS
        # for i in range(self.columns_out.shape[0]):
            # F = F - LoadFactor*AppliedDirichlet[i]*stiffness.getcol(self.columns_out[i])

        # MUCH FASTER APPROACH
        # F = F - (stiffness[:,self.columns_out]*AppliedDirichlet*LoadFactor)[:,None]
        nnz_cols = ~np.isclose(AppliedDirichlet,0.0)
        if self.columns_out[nnz_cols].shape[0]==0:
            F[self.columns_in] = F[self.columns_in]
        else:
            F[self.columns_in] = F[self.columns_in] - (stiffness[self.columns_in,:]\
                [:,self.columns_out[nnz_cols]]*AppliedDirichlet[nnz_cols]*LoadFactor)[:,None]

        if only_residual:
            return F

        # GET REDUCED FORCE VECTOR
        F_b = F[self.columns_in,0]

        # GET REDUCED STIFFNESS
        stiffness_b = stiffness[self.columns_in,:][:,self.columns_in]

        # GET REDUCED MASS MATRIX
        if self.analysis_type != 'static':
            mass_b = mass[self.columns_in,:][:,self.columns_in]
            return stiffness_b, F_b, F, mass_b

        return stiffness_b, F_b, F

    def UpdateFixDoFs(self, AppliedDirichletInc, fsize, nvar):
        """Updates the geometry (DoFs) with incremental Dirichlet boundary conditions
            for fixed/constrained degrees of freedom only. Needs to be applied per time steps"""

        # GET TOTAL SOLUTION
        TotalSol = np.zeros((fsize,1))
        TotalSol[self.columns_out,0] = AppliedDirichletInc

        # RE-ORDER SOLUTION COMPONENTS
        dU = TotalSol.reshape(int(TotalSol.shape[0]/nvar),nvar)

        return dU

    def UpdateFreeDoFs(self, sol, fsize, nvar):
        """Updates the geometry with iterative solutions of Newton-Raphson
            for free degrees of freedom only. Needs to be applied per time NR iteration"""

        # GET TOTAL SOLUTION
        TotalSol = np.zeros((fsize,1))
        TotalSol[self.columns_in,0] = sol

        # RE-ORDER SOLUTION COMPONENTS
        dU = TotalSol.reshape(int(TotalSol.shape[0]/nvar),nvar)

        return dU
