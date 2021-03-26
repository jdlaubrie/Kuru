import numpy as np
from Kuru import QuadratureRule, FunctionSpace , Mesh
from Kuru.FiniteElements.LocalAssembly._KinematicMeasures_ import _KinematicMeasures_
from Kuru.VariationalPrinciple._GeometricStiffness_ import GeometricStiffnessIntegrand as GetGeomStiffness
from .DisplacementApproachIndices import FillGeometricB
#from ._MassIntegrand_ import __MassIntegrand__, __ConstantMassIntegrand__

__all__ = ["VariationalPrinciple"]

class VariationalPrinciple(object):

    energy_dissipation = []
    internal_energy = []
    kinetic_energy = []
    external_energy = []

    power_dissipation = []
    internal_power = []
    kinetic_power = []
    external_power = []


    def __init__(self, mesh, variables_order=(1,0),
        analysis_type='static', analysis_nature='nonlinear', fields='mechanics',
        quadrature_rules=None, median=None, quadrature_type=None,
        function_spaces=None, compute_post_quadrature=True):

        self.variables_order = variables_order
        self.nvar = None
        self.ndim = mesh.points.shape[1]

        if isinstance(self.variables_order,int):
            self.variables_order = tuple(self.variables_order)

        self.quadrature_rules = quadrature_rules
        self.quadrature_type = quadrature_type
        self.function_spaces = function_spaces
        self.median = median
        self.analysis_type = analysis_type
        self.analysis_nature = analysis_nature
        self.fields = fields

        self.compute_post_quadrature = compute_post_quadrature

        # GET NUMBER OF VARIABLES
        self.GetNumberOfVariables()

    def GetQuadratureOrder(self, C, element_type, quadrature_degree=None):
        """Finds quadrature degree/strength for a given polynomial order C=p-1 [where p is polynomial degree]"""
        if quadrature_degree is None:
            if element_type == "tri" or element_type == "tet":
                norder = 2*C if C > 0 else 1
                norder_post = 2*(C+1)
            else:
                norder = C+2
                # ACTUAL
                # norder_post = 2*(C+2)
                # ALTHOUGH THIS INTEGRATES EXACTLY
                norder_post = C+2
        else:
            norder = quadrature_degree
            if element_type == "tri" or element_type == "tet":
                norder_post = 2*quadrature_degree
            else:
                norder_post = quadrature_degree

        return norder, norder_post

    def GetQuadraturesAndFunctionSpaces(self, mesh, variables_order=(1,),
        quadrature_rules=None, quadrature_type=None, function_spaces=None, compute_post_quadrature=True,
        equally_spaced_bases=False, quadrature_degree=None):
        """"The default function for computing quadrature rules and function spaces for equall order single
            and multi-physics/fields problems"""

        C = mesh.InferPolynomialDegree() - 1
        mesh.InferBoundaryElementType()

        if quadrature_rules == None and self.quadrature_rules == None:

            # OPTION FOR QUADRATURE TECHNIQUE FOR TRIS AND TETS
            optimal_quadrature = 3
            if mesh.element_type == "quad" or mesh.element_type == "hex":
                if quadrature_type == "wv":
                    optimal_quadrature = 4

            norder, norder_post = self.GetQuadratureOrder(C, mesh.element_type, quadrature_degree=quadrature_degree)

            # GET QUADRATURE
            quadrature = QuadratureRule(optimal=optimal_quadrature, norder=norder, mesh_type=mesh.element_type)
            if self.compute_post_quadrature:
                # COMPUTE INTERPOLATION FUNCTIONS AT ALL INTEGRATION POINTS FOR POST-PROCESSING
                post_quadrature = QuadratureRule(optimal=optimal_quadrature, norder=norder_post, mesh_type=mesh.element_type)
            else:
                post_quadrature = None

            # BOUNDARY QUADRATURE
            bquadrature = QuadratureRule(optimal=optimal_quadrature, norder=C+2, mesh_type=mesh.boundary_element_type)

            self.quadrature_rules = (quadrature,post_quadrature,bquadrature)
        else:
            self.quadrature_rules = quadrature_rules

        if function_spaces == None and self.function_spaces == None:

            # CREATE FUNCTIONAL SPACES
            function_space = FunctionSpace(mesh, self.quadrature_rules[0], p=C+1, equally_spaced=equally_spaced_bases)
            if self.compute_post_quadrature:
                post_function_space = FunctionSpace(mesh, self.quadrature_rules[1], p=C+1, equally_spaced=equally_spaced_bases)
            else:
                post_function_space = None

            # CREATE BOUNDARY FUNCTIONAL SPACES
            bfunction_space = FunctionSpace(mesh.CreateDummyLowerDimensionalMesh(),
                self.quadrature_rules[2], p=C+1, equally_spaced=equally_spaced_bases)

            self.function_spaces = (function_space,post_function_space,bfunction_space)
        else:
            self.function_spaces = function_spaces


        local_size = self.function_spaces[0].Bases.shape[0]*self.nvar
        self.local_rows = np.repeat(np.arange(0,local_size),local_size,axis=0)
        self.local_columns = np.tile(np.arange(0,local_size),local_size)
        self.local_size = local_size

        # FOR MASS
        local_size_m = self.function_spaces[0].Bases.shape[0]*self.ndim
        self.local_rows_mass = np.repeat(np.arange(0,local_size_m),local_size_m,axis=0)
        self.local_columns_mass = np.tile(np.arange(0,local_size_m),local_size_m)
        self.local_size_m = local_size_m

    def GetNumberOfVariables(self):
        """Returns (self.nvar) i.e. number of variables/unknowns per node, for the formulation.
            Note that self.nvar does not take into account the unknowns which get condensated
        """

        # nvar = 0
        # for i in self.variables_order:
        #     # DO NOT COUNT VARIABLES THAT GET CONDENSED OUT
        #     if i!=0:
        #         if mesh.element_type == "tri":
        #             nvar += (i+1)*(i+2) // 2
        #         elif mesh.element_type == "tet":
        #             nvar += (i+1)*(i+2)*(i+3) // 6
        #         elif mesh.element_type == "quad":
        #             nvar += (i+1)**2
        #         elif mesh.element_type == "hex":
        #             nvar += (i+1)**3

        # nvar = sum(self.variables_order)
        if self.nvar == None:
            self.nvar = self.ndim
        return self.nvar

    def FindIndices(self,A):
        return self.local_rows, self.local_columns, A.ravel()

    def GeometricStiffnessIntegrand(self, SpatialGradient, CauchyStressTensor):
        """Applies to displacement based, displacement potential based and all mixed
        formulations that involve static condensation"""

        ndim = self.ndim
        nvar = self.nvar

        B = np.zeros((nvar*SpatialGradient.shape[0],ndim*ndim))
        S = np.zeros((ndim*ndim,ndim*ndim))
        SpatialGradient = SpatialGradient.T.copy('c')

        FillGeometricB(B,SpatialGradient,S,CauchyStressTensor,ndim,nvar)

        BDB = np.dot(np.dot(B,S),B.T)

        return BDB


    def __GeometricStiffnessIntegrand__(self, SpatialGradient, CauchyStressTensor, detJ):
        """Applies to displacement based formulation"""
        return GetGeomStiffness(np.ascontiguousarray(SpatialGradient),CauchyStressTensor, detJ, self.nvar)

    def VolumetricStiffnessIntegrand(self, material, SpatialGradient, detJ, dV):
        """Computes the volumetric stiffness using Hu-Washizu on Mean Dilatation method"""

        if material.has_low_level_dispatcher:
            from ._VolumetricStiffness_ import _VolumetricStiffnessIntegrand_
            stiffness, MeanVolume = _VolumetricStiffnessIntegrand_(material,
                np.ascontiguousarray(SpatialGradient), np.ascontiguousarray(detJ),
                np.ascontiguousarray(dV), self.nvar)
        else:
            MaterialVolume = np.sum(dV)
            if material.has_state_variables and material.has_growth_remodeling:
                dve = np.true_divide(detJ,material.StateVariables[:,material.id_growth])
                CurrentElasticVolume = np.sum(dve)
                # AVERAGE SPATIAL GRADIENT IN PHYSICAL ELEMENT [\frac{1}{v}\int\nabla(N)dv(nodeperelem x ndim)]
                AverageDeformationv = np.einsum('i,ijk,i->jk',material.StateVariables[:,material.id_density],SpatialGradient,dve)
                AverageDeformationv = AverageDeformationv.flatten()
                AverageDeformationu = np.einsum('ijk,i->jk',SpatialGradient,dve)
                AverageDeformationu = AverageDeformationu.flatten()
                stiffness = np.einsum('i,j->ij',AverageDeformationv,AverageDeformationu)
                MeanVolume = (CurrentElasticVolume-MaterialVolume)/MaterialVolume
            elif material.has_state_variables and not material.has_growth_remodeling:
                CurrentElasticVolume = np.sum(detJ)
                # AVERAGE SPATIAL GRADIENT IN PHYSICAL ELEMENT [\frac{1}{v}\int\nabla(N)dv(nodeperelem x ndim)]
                AverageDeformationv = np.einsum('i,ijk,i->jk',material.StateVariables[:,material.id_density],SpatialGradient,detJ)
                AverageDeformationv = AverageDeformationv.flatten()
                AverageDeformationu = np.einsum('ijk,i->jk',SpatialGradient,detJ)
                AverageDeformationu = AverageDeformationu.flatten()
                stiffness = np.einsum('i,j->ij',AverageDeformationv,AverageDeformationu)
                MeanVolume = (CurrentElasticVolume-MaterialVolume)/MaterialVolume
            elif not material.has_state_variables and not material.has_growth_remodeling:
                CurrentVolume = np.sum(detJ)
                # AVERAGE SPATIAL GRADIENT IN PHYSICAL ELEMENT [\frac{1}{v}\int\nabla(N)dv(nodeperelem x ndim)]
                AverageSpatialGradient = np.einsum('ijk,i->jk',SpatialGradient,detJ)
                AverageSpatialGradient = AverageSpatialGradient.flatten()
                stiffness = np.einsum('i,j->ij',AverageSpatialGradient,AverageSpatialGradient)
                MeanVolume = (CurrentVolume-MaterialVolume)/MaterialVolume

            stiffness = np.true_divide(stiffness,MaterialVolume)

        material.pressure = material.kappa*MeanVolume
        stiffness *= material.kappa

        return stiffness



