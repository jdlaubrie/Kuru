# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Kuru
sys.path.append("/media/jdlaubrie/c57cebce-8099-4a2e-9731-0bf5f6e163a5/Kuru/")
#sys.path.append(os.path.expanduser("~"))
from Kuru import *

#============================================================
#============= MATERIALS SETS  ================
#============================================================
def GetElementSets(mesh,no_materials):
    """ Organize elements by set, then each of those set will be link
        with a material"""

    print("Getting elements and nodes set for materials ...")

    # ELEMENTS AND NODE SETS
    # arterial wall {1,2}
    element_sets = [[] for i in range(no_materials)]
    element_sets[0] = np.hstack((np.where(mesh.element_to_set==0)[0],np.where(mesh.element_to_set==1)[0]))

    node_sets = [[] for i in range(no_materials)]
    for elset in range(no_materials):
        node_sets[elset] = np.unique(mesh.elements[element_sets[elset]])

    return element_sets, node_sets

#============================================================
#============= ANISOTROPIC FIBRE DIRECTIONS  ================
#============================================================
def Directions(mesh):
    """
        Routine dedicated to compute the fibre direction of components by node for 
        the anisotropic Material in Florence. First three directions are taken into
        the code for Rotation matrix, so always it should be present in this order,
        Normal, Tangential, Axial.
    """

    print("Building anisotrpic directions ...")

    ndim = mesh.InferSpatialDimension()
    nfibre = 6
    # Geometric definitions per element
    divider = mesh.elements.shape[1]
    directrix = [0.,0.,1.]
    direction = np.zeros((mesh.nelem,nfibre,ndim),dtype=np.float64)
    # Loop throught the element in the mesh
    for elem in range(mesh.nelem):
        # Geometric definitions per element
        center = np.sum(mesh.points[mesh.elements[elem,:],:],axis=0)/divider
        tangential = np.cross(directrix,center)
        tangential = tangential/np.linalg.norm(tangential)
        normal = np.cross(tangential,directrix)
        direction[elem][0][:] = normal/np.linalg.norm(normal)
        # Define the anisotropic orientations
        direction[elem][1][:]=tangential
        direction[elem][2][:]=directrix
        direction[elem][3][:]=np.multiply(directrix,np.cos(np.pi/4.)) + np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][4][:]=np.multiply(directrix,np.cos(np.pi/4.)) - np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][5][:]=tangential

    return direction

#============================================================
#===============  HOMOGENIZED CMT  ==========================
#============================================================
def CylinderDepositionStretch(p=1,growth=False):

    # build the path to the mesh file
    ProblemPath = os.getcwd()
    mesh_file = ProblemPath + '/Cylinder.msh'

    #===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    ndim = mesh.InferSpatialDimension()
    mesh.GetHighOrderMesh(p=p)

    #Boolean arrays for boundary condition in Dirichlet
    BottomSurface = np.zeros(mesh.nnode,dtype=bool)
    TopSurface = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_Y = np.zeros(mesh.nnode,dtype=bool)
    InnerSurface = np.zeros(mesh.faces.shape[0],dtype=bool)
    OuterSurface = np.zeros(mesh.faces.shape[0],dtype=bool)

    # bottom surface, Z=0. {5};
    BottomSurface[mesh.faces[np.where(mesh.face_to_surface==4)[0],:]] = True
    # top surface, Z=length. {49};
    TopSurface[mesh.faces[np.where(mesh.face_to_surface==48)[0],:]] = True
    # symmetry X surface. {22,44};
    Symmetry_X[mesh.faces[np.where(mesh.face_to_surface==21)[0],:]] = True
    Symmetry_X[mesh.faces[np.where(mesh.face_to_surface==43)[0],:]] = True
    # symmetry Y surface. {14,36};
    Symmetry_Y[mesh.faces[np.where(mesh.face_to_surface==13)[0],:]] = True
    Symmetry_Y[mesh.faces[np.where(mesh.face_to_surface==35)[0],:]] = True
    # inner surface. {26,48};
    InnerSurface[np.where(mesh.face_to_surface==25)[0]] = True
    InnerSurface[np.where(mesh.face_to_surface==47)[0]] = True
    # outer surface. {18,40};
    OuterSurface[np.where(mesh.face_to_surface==17)[0]] = True
    OuterSurface[np.where(mesh.face_to_surface==39)[0]] = True

    DirichletBoundary = {}
    DirichletBoundary['Bottom'] = BottomSurface
    DirichletBoundary['Top'] = TopSurface
    DirichletBoundary['SymmetryX'] = Symmetry_X
    DirichletBoundary['SymmetryY'] = Symmetry_Y
    RobinBoundary = {}
    RobinBoundary['Inner'] = InnerSurface
    RobinBoundary['Outer'] = OuterSurface

    #===============  MATERIAL DEFINITION  ====================
    no_materials = 1
    # Get sets for each material
    element_sets, node_sets = GetElementSets(mesh,no_materials)

    # Set the state variables for each set
    state_variables = np.zeros((node_sets[0].shape[0],21),dtype=np.float64)
    # Remodeling is the inverse of Deposition Stretches
    state_variables[:,0] = 1.34*1.25
    state_variables[:,4] = 1./1.34
    state_variables[:,8] = 1./1.25
    state_variables[:,9] = 1./1.1
    state_variables[:,10] = 1./1.062
    state_variables[:,11] = 1./1.062
    state_variables[:,12] = 1./1.062
    state_variables[:,13] = 1./1.062
    # Total initial density
    state_variables[:,14] = 241.5
    state_variables[:,15] = 157.5
    state_variables[:,16] = 65.1
    state_variables[:,17] = 260.4
    state_variables[:,18] = 260.4
    state_variables[:,19] = 65.1
    state_variables[:,20] = 1.0

    # fibre directions [thick,smc,co1,co2,co3,co4]
    fibre_direction = Directions(mesh)
    
    # material elastic properties
    mu = 72.0e-6
    k1m = 7.6e-6
    k2m = 11.4
    k1c = 568.0e-6
    k2c = 11.2

    # Define hyperelastic material for media
    artery = ArterialWallMixture(ndim,
            is_nearly_incompressible=False,
            id_growth=20,
            id_density=14,
            rho0=1050.,
            mu=mu,
            kappa=mu*100.0,
            k1m=k1m,
            k2m=k2m,
            k1c=k1c,
            k2c=k2c,
            maxi_active_stress=54.0e-3,
            maxi_active_stretch=1.4,
            zero_active_stretch=0.8,
            active_stretch=1.0,
            anisotropic_orientations=fibre_direction,
            #load_factor=[0.,0.,1.],
            state_variables=state_variables,
            element_set=element_sets[0],
            node_set=node_sets[0])

    # Material tuple
    materials = (artery, )
    # kappa/mu=20  => nu=0.475 (Poisson's ratio)
    # kappa/mu=33  => nu=0.485 (Poisson's ratio)
    # kappa/mu=100 => nu=0.495 (Poisson's ratio)

    #==================  FORMULATION  =========================
    formulation = DisplacementFormulation(mesh)

    #===============  BOUNDARY CONDITIONS  ====================
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, DirichletBoundary):
        boundary_data = np.zeros((mesh.nnode, 3))+np.NAN
        # boundary conditions base on BoundarySurface boolean array
        boundary_data[DirichletBoundary['Bottom'],2] = 0.
        boundary_data[DirichletBoundary['Top'],2] = 0.
        boundary_data[DirichletBoundary['SymmetryX'],0] = 0.
        boundary_data[DirichletBoundary['SymmetryY'],1] = 0.

        return boundary_data

    # Pressure Boundary Conditions
    def Robin_Function(mesh, RobinBoundary):
        # Dicts = {type: string, flags: logical_array, data: float_array}
        # the type would be Pressure, Spring or Dashpot
        Pressure = {'type': 'Pressure'}
        Pressure['flags'] = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        Pressure['data'] = np.zeros((mesh.faces.shape[0]))
        # Inner pressure
        Pressure['flags'][RobinBoundary['Inner']] = True
        Pressure['data'][RobinBoundary['Inner']] = -13.3322e-3

        #Spring = {'type': 'Spring'}
        #Spring['flags'] = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        #Spring['data'] = np.zeros((mesh.faces.shape[0]))
        # Outer spring
        #Spring['flags'][RobinBoundary['Outer']] = True
        #Spring['data'][RobinBoundary['Outer']] = 10.0e-6

        return Pressure #, Spring

    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetRobinCriteria(Robin_Function, mesh, RobinBoundary)

    #================= SOLUTION  =======================
    if not growth:
        # set solver parameters
        fem_solver = FEMSolver(analysis_nature="nonlinear",
                           analysis_type="static",
                           break_at_stagnation=False,
                           maximum_iteration_for_newton_raphson=50,
                           optimise=False,
                           parallelise=False,
                           print_incremental_log=True,
                           has_moving_boundary=True,
                           #load_factor=[0.34,0.33,0.33],
                           number_of_load_increments=1)

        # Call the solver
        solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
            materials=materials, boundary_condition=boundary_condition)

    else:
        # set solver parameters for growth and remodeling. 60 days
        fem_solver = FEMSolver(analysis_nature="nonlinear",
                           analysis_type="static",
                           break_at_stagnation=False,
                           maximum_iteration_for_newton_raphson=50,
                           optimise=False,
                           parallelise=False,
                           print_incremental_log=True,
                           has_moving_boundary=True,
                           #load_factor=[0.34,0.33,0.33],
                           number_of_time_increments=6,
                           total_time=60.)

        # growth and remodeling time integratos
        growth_remodeling = ExplicitGrowthRemodelingIntegrator(gain=0.05,
                               turnover=101.0,
                               density_turnover="self",
                               damage_spread_space=10.0,
                               damage_axis=2)

        # Call the solver
        solution = fem_solver.Solve(formulation=formulation, mesh=mesh, materials=materials,
            boundary_condition=boundary_condition, growth_remodeling=growth_remodeling)


    # Write VTK file for visualization
    solution.WriteVTK('Cylinder', quantity=0, time_problem=False) #disp #42-s_VM

if __name__ == "__main__":
    CylinderDepositionStretch(p=1,growth=True)
    
