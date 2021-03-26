import numpy as np
from .LineBP import LagrangeBP_, LagrangeGaussLobattoBP_

def Lagrange(C,xi):

    n = C+2
    ranger = np.arange(n)
    eps = np.linspace(-1.,1.,n)

    A = np.zeros((n,n))
    A[:,0] = 1.
    for i in range(1,n):
        A[:,i] = eps**i

    RHS = np.zeros((n,n))
    np.fill_diagonal(RHS,1)
    coeff = np.linalg.solve(A,RHS)
    xis = np.ones(n)*xi**ranger
    N = np.dot(coeff.T,xis)
    dN = np.dot(coeff[1:,:].T,xis[:-1]*(1+ranger[:-1]))

    return (N,dN,eps)

def LagrangeGaussLobatto(C,xi):
    
    from Kuru.QuadratureRules import GaussLobattoQuadrature

    n = C+2
    ranger = np.arange(n)

    eps = GaussLobattoQuadrature(n)[0][:,0]

    A = np.zeros((n,n))
    A[:,0] = 1.
    for i in range(1,n):
        A[:,i] = eps**i
    # A1[:,1:] = np.array([eps**i for i in range(1,n)]).T[0,:,:]


    RHS = np.zeros((n,n))
    np.fill_diagonal(RHS,1)
    coeff = np.linalg.solve(A,RHS)
    xis = np.ones(n)*xi**ranger
    N = np.dot(coeff.T,xis)
    # dN = np.einsum('i,ij,i',1+ranger[:-1],coeff[1:,:],xis[:-1])
    dN = np.dot(coeff[1:,:].T,xis[:-1]*(1+ranger[:-1]))

    return (N,dN,eps)

def Legendre(C,xi):
    # For Linear Basis Generating Legendre Polynomials is Not Required
    if C==0:
        N = np.array([1.0/2*(1-xi), 1.0/2*(1+xi)])
        dN = np.array([-1.0/2, 1./2])

    # For Higher Order 
    elif C>0:
        # The first two Legendre polynomials 
        p0 = 1.0; p1 = xi
        # Derivatives of The First Two Legendre Polynomials 
        dp0 = 0.0; dp1 = 1.0
        # Allocate size and dimensions
        ndim = C+2
        P = np.zeros((ndim+1,1)); dP = np.zeros((ndim+1,1))
        N = np.zeros((ndim+1,1)); dN = np.zeros((ndim+1,1))
        P[0] = p0; P[1] = p1
        dP[0] = dp0; dP[1] = dp1
        # Generate Legendre Polynomials
        for i in range(2,ndim+1):
            P[i]  = ((2.0*i-1)*xi*P[i-1] - (i-1)*P[i-2])/(i)
            dP[i]  = ((2.0*i-1)*xi*dP[i-1] + (2.0*i-1)*P[i-1] - (i-1)*dP[i-2])/(i)

        # From Legendre Polynomials Generate FE Basis Functions 
        for i in range(3,ndim+2):
            # N[i-1] =  (P[i-1]-P[i-3])/np.sqrt(2*(2*i-3))
            # dN[i-1] =  (dP[i-1]-dP[i-3])/np.sqrt(2*(2*i-3))
            # Ledger's Normalisation 
            N[i-1] =  (P[i-1]-P[i-3])/((2.0*i-3.))
            dN[i-1] =  (dP[i-1]-dP[i-3])/((2.0*i-3.))


        # Put the hat functions at exterior nodes  
        N = np.append([np.append([1.0/2.0*(1.0-xi)],N[2:-1])],[1.0/2*(1.0+xi)])
        dN = np.append([np.append([-0.5],dN[2:-1])],[0.5])


    return (N,dN)

def LagrangeBP(C,xi):
    return LagrangeBP_(C,xi)

def LagrangeGaussLobattoBP(C,xi):
    
    n = C+2
    from Kuru.QuadratureRules import GaussLobattoQuadrature
    eps = GaussLobattoQuadrature(n)[0][:,0].copy()
    return LagrangeGaussLobattoBP_(C,xi,eps)
