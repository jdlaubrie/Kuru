#include "_MaterialBase_.h"

template<typename U>
class _NeoHookean_2_ : public _MaterialBase_<U> {
public:
    U mu;
    U kappa;

    FASTOR_INLINE _NeoHookean_2_() = default;

    FASTOR_INLINE
    _NeoHookean_2_(U mu, U kappa) {
        this->mu = mu;
        this->kappa = kappa;
    }

    FASTOR_INLINE
    void SetParameters(U mu, U kappa){
        this->mu = mu;
        this->kappa = kappa;
    }


    template<typename T=U, size_t ndim>
    FASTOR_INLINE
    std::tuple<Tensor<T,ndim,ndim>, typename MechanicsHessianType<T,ndim>::return_type>
    _KineticMeasures_(const T *Fnp, int near_incomp, const T pressure) {

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
        if (ndim==2) {trb += 1.;}
        T coeff0 = std::pow(J,-5./3.);

        // ISOCHORIC PART
        Tensor<T,ndim,ndim> stress = mu*coeff0*(b-1./3.*trb*I);

        // FIND ELASTICITY TENSOR
        auto II_ijkl = einsum<Index<i,j>,Index<k,l>>(I,I);
        auto II_ikjl = permutation<Index<i,k,j,l>>(II_ijkl);
        auto II_iljk = permutation<Index<i,l,j,k>>(II_ijkl);

        auto Ib_ijkl = einsum<Index<i,j>,Index<k,l>>(I,b);
        auto bI_ijkl = einsum<Index<i,j>,Index<k,l>>(b,I);

        // ISOCHORIC PART
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

        auto hessian = voigt(elasticity);

        auto kinetics = std::make_tuple(stress,hessian);
        return kinetics;
    }



    template<typename T>
    void KineticMeasures(T *Snp, T* Hnp, int ndim, int ngauss, const T *Fnp, int near_incomp, const T pressure);

};

template<> template<>
void _NeoHookean_2_<Real>::KineticMeasures<Real>(Real *Snp, Real* Hnp,
    int ndim, int ngauss, const Real *Fnp, int near_incomp, const Real pressure) {

    if (ndim==3) {
        Tensor<Real,3,3> stress;
        Tensor<Real,6,6> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,3>(Fnp+9*g, near_incomp, pressure);
            copy_fastor(Snp,stress,g*9);
            copy_fastor(Hnp,hessian,g*36);
        }
    }
    else if (ndim==2) {
        Tensor<Real,2,2> stress;
        Tensor<Real,3,3> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,2>(Fnp+4*g, near_incomp, pressure);
            copy_fastor(Snp,stress,g*4);
            copy_fastor(Hnp,hessian,g*9);
        }
    }
}
