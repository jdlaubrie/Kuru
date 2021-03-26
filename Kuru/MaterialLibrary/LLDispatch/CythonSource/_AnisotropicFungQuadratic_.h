#include "_MaterialBase_.h"

template<typename U>
class _AnisotropicFungQuadratic_ : public _MaterialBase_<U> {
public:
    U mu;
    U kappa;
    U k1;
    U k2;

    FASTOR_INLINE _AnisotropicFungQuadratic_() = default;

    FASTOR_INLINE
    _AnisotropicFungQuadratic_(U mu, U kappa, U k1, U k2) {
        this->mu = mu;
        this->kappa = kappa;
        this->k1 = k1;
        this->k2 = k2;
    }

    FASTOR_INLINE
    void SetParameters(U mu, U kappa, U k1, U k2){
        this->mu = mu;
        this->kappa = kappa;
        this->k1 = k1;
        this->k2 = k2;
    }

    template<typename T=U, size_t ndim>
    FASTOR_INLINE
    std::tuple<Tensor<T,ndim,ndim>, typename MechanicsHessianType<T,ndim>::return_type>
    _KineticMeasures_(const T *Fnp, const T *Nnp, int nfibre, int near_incomp, const T pressure) {

        // CREATE FASTOR TENSORS
        Tensor<T,ndim,ndim> F;
        // COPY NUMPY ARRAY TO FASTOR TENSORS
        copy_numpy(F,Fnp);

        // FIND THE KINEMATIC MEASURES
        Tensor<Real,ndim,ndim> I; I.eye2();
        auto J = determinant(F);
        auto b = matmul(F,transpose(F));

        // COMPUTE CAUCHY STRESS TENSOR
        T trb = trace(b);
        if (ndim==2) { trb += 1.; }
        T coeff0 = std::pow(J,-5./3.);

        // ISOCHORIC PART OF STRESS
        Tensor<T,ndim,ndim> stress = mu*coeff0*(b-1./3.*trb*I);

        // FIND ELASTICITY TENSOR
        auto II_ijkl = einsum<Index<i,j>,Index<k,l>>(I,I);
        auto II_ikjl = permutation<Index<i,k,j,l>>(II_ijkl);
        auto II_iljk = permutation<Index<i,l,j,k>>(II_ijkl);

        auto Ib_ijkl = einsum<Index<i,j>,Index<k,l>>(I,b);
        auto bI_ijkl = einsum<Index<i,j>,Index<k,l>>(b,I);

        // ISOCHORIC PART OF ELASTICITY TENSOR
        Tensor<T,ndim,ndim,ndim,ndim> elasticity = 2.*mu*coeff0*(1./9.*trb*II_ijkl-
            1./3.*(Ib_ijkl+bI_ijkl)+1./6.*trb*(II_ikjl+II_iljk));

        // VOLUMETRIC PART
        if (near_incomp) {
            stress += pressure*I;
            elasticity += pressure*(II_ijkl-(II_ikjl+II_iljk));
        }
        else {
            stress += kappa*(J-1.)*I;
            elasticity += kappa*((2.*J-1.)*II_ijkl-(J-1.)*(II_ikjl+II_iljk));
        }

        // LOOP OVER FIBRES ORIENTATIONS
        for (int n=0; n<nfibre; ++n) {
            Tensor<T,ndim> N;
            copy_numpy(N,Nnp+3*n);
            auto FN = matmul(F,N);
            auto innerFN = inner(FN,FN);
            auto outerFN = outer(FN,FN);
            T coeff1 = std::pow((innerFN-1.),2);
            T k1c = ((innerFN-1.0)>=0.0) ? k1 : 0.075*k1;
            T coeff2 = std::exp(k2*coeff1);
            // Cauchy stress of fibre
            stress += 2.*k1c/J*(innerFN-1.)*coeff2*outerFN;
            // Elasticity tensor of fibre
            auto FNFN_ijkl = einsum<Index<i,j>,Index<k,l>>(outerFN,outerFN);
            elasticity += 4.*k1c/J*(1.+2.*k2*coeff1)*coeff2*FNFN_ijkl;
        }

        auto hessian = voigt(elasticity);

        auto kinetics = std::make_tuple(stress,hessian);
        return kinetics;
    }



    template<typename T>
    void KineticMeasures(T *Snp, T* Hnp, int ndim, int ngauss, const T *Fnp, int nfibre, const T *Nnp, int near_incomp, const T pressure);

};

template<> template<>
void _AnisotropicFungQuadratic_<Real>::KineticMeasures<Real>(Real *Snp, Real* Hnp,
    int ndim, int ngauss, const Real *Fnp, int nfibre, const Real *Nnp, int near_incomp, const Real pressure) {

    if (ndim==3) {
        Tensor<Real,3,3> stress;
        Tensor<Real,6,6> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,3>(Fnp+9*g, Nnp, nfibre, near_incomp, pressure);
            copy_fastor(Snp,stress,g*9);
            copy_fastor(Hnp,hessian,g*36);
        }
    }
    else if (ndim==2) {
        Tensor<Real,2,2> stress;
        Tensor<Real,3,3> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,2>(Fnp+4*g, Nnp, nfibre, near_incomp, pressure);
            copy_fastor(Snp,stress,g*4);
            copy_fastor(Hnp,hessian,g*9);
        }
    }
}
