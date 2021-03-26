#ifndef _HELPER__H
#define _HELPER__H

#include <Fastor/Fastor.h>
#include <tuple>
using namespace Fastor;

typedef double Real;

enum 
{
    i,j,k,l,m,n,o,p,q,r
};

template<typename T, size_t ... Dims> 
FASTOR_INLINE
void copy_numpy(Tensor<T,Dims...> &A, const T* A_np, size_t offset=0) {
    std::copy(A_np,A_np+A.size(),A.data());
}


template<typename T, size_t ... Dims> 
FASTOR_INLINE
void copy_fastor(T* A_np, const Tensor<T,Dims...> &A, size_t offset=0) {
    std::copy(A.data(),A.data()+A.size(),A_np+offset);
}


// For mechanics
template<typename T, size_t N>
struct MechanicsHessianType;
template<typename T>
struct MechanicsHessianType<T,2> {
    using return_type = Tensor<T,3,3>;
};
template<typename T>
struct MechanicsHessianType<T,3> {
    using return_type = Tensor<T,6,6>;
};

// For possion
template<typename T, size_t N>
struct PoissonHessianType;
template<typename T>
struct PoissonHessianType<T,2> {
    using return_type = Tensor<T,2,2>;
};
template<typename T>
struct PoissonHessianType<T,3> {
    using return_type = Tensor<T,3,3>;
};

#endif
