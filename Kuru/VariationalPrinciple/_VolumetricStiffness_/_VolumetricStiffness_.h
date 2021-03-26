#ifndef VOLUMETRICSTIFFNESS_H
#define VOLUMETRICSTIFFNESS_H

#include <algorithm>
#include <cstdint>
#include <cstdlib>

#ifndef LL_TYPES
#define LL_TYPES
using Real = double;
using Integer = std::int64_t;
using UInteger = std::uint64_t;
#endif

void _VolumetricStiffnessFiller_(Real *volumetric_stiffness, Real *mean_volume,
    const Real *SpatialGradient, const Real *detJ, const Real *dV,
    const Real *density, const Real *Growth, int has_stvar, int has_gr,
    Integer ndim, Integer nvar, Integer nodeperelem, Integer ngauss)
{
    if (has_stvar && has_gr) {

    Real *dve   = (Real*)malloc(sizeof(Real)*ngauss);
    Real volumex = 0.;
    Real VolumeX = 0.;
    for (Integer g=0; g<ngauss; ++g) {
        dve[g] = detJ[g]/Growth[g];
        volumex += dve[g];
        VolumeX += dV[g];
    }
    mean_volume[0] = (volumex-VolumeX)/VolumeX;

    Real *AverageSpatialGradientv   = (Real*)malloc(sizeof(Real)*ndim*nodeperelem);
    Real *AverageSpatialGradientu   = (Real*)malloc(sizeof(Real)*ndim*nodeperelem);
    std::fill(AverageSpatialGradientv,AverageSpatialGradientv+nodeperelem*ndim,0.);
    std::fill(AverageSpatialGradientu,AverageSpatialGradientu+nodeperelem*ndim,0.);
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer i=0; i<ndim; ++i) {
            for (Integer g=0; g<ngauss; ++g) {
              AverageSpatialGradientv[a*ndim+i] += density[g]*SpatialGradient[g*nodeperelem*ndim+a*ndim+i]*dve[g];
              AverageSpatialGradientu[a*ndim+i] += SpatialGradient[g*nodeperelem*ndim+a*ndim+i]*dve[g];
            }
        }
    }
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer b=0; b<nodeperelem; ++b) {
            for (Integer i=0; i<ndim; ++i) {
                for (Integer j=0; j<ndim; ++j) {
                    volumetric_stiffness[(a*nvar+i)*nodeperelem*nvar+(b*nvar+j)] += 
                             AverageSpatialGradientv[a*ndim+i]*AverageSpatialGradientu[b*ndim+j]/VolumeX;
                }
            }
        }
    }
    free(dve);
    free(AverageSpatialGradientu);
    free(AverageSpatialGradientv);
    }
    
    else if (has_stvar && !(has_gr)) {

    Real volumex = 0.;
    Real VolumeX = 0.;
    for (Integer g=0; g<ngauss; ++g) {
        volumex += detJ[g];
        VolumeX += dV[g];
    }
    mean_volume[0] = (volumex-VolumeX)/VolumeX;

    Real *AverageSpatialGradientv   = (Real*)malloc(sizeof(Real)*ndim*nodeperelem);
    Real *AverageSpatialGradientu   = (Real*)malloc(sizeof(Real)*ndim*nodeperelem);
    std::fill(AverageSpatialGradientv,AverageSpatialGradientv+nodeperelem*ndim,0.);
    std::fill(AverageSpatialGradientu,AverageSpatialGradientu+nodeperelem*ndim,0.);
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer i=0; i<ndim; ++i) {
            for (Integer g=0; g<ngauss; ++g) {
              AverageSpatialGradientv[a*ndim+i] += density[g]*SpatialGradient[g*nodeperelem*ndim+a*ndim+i]*detJ[g];
              AverageSpatialGradientu[a*ndim+i] += SpatialGradient[g*nodeperelem*ndim+a*ndim+i]*detJ[g];
            }
        }
    }
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer b=0; b<nodeperelem; ++b) {
            for (Integer i=0; i<ndim; ++i) {
                for (Integer j=0; j<ndim; ++j) {
                    volumetric_stiffness[(a*nvar+i)*nodeperelem*nvar+(b*nvar+j)] += 
                             AverageSpatialGradientv[a*ndim+i]*AverageSpatialGradientu[b*ndim+j]/VolumeX;
                }
            }
        }
    }
    free(AverageSpatialGradientu);
    free(AverageSpatialGradientv);
    }

    else if (!(has_stvar) && !(has_gr)) {

    Real volumex = 0.;
    Real VolumeX = 0.;
    for (Integer g=0; g<ngauss; ++g) {
        volumex += detJ[g];
        VolumeX += dV[g];
    }
    mean_volume[0] = (volumex-VolumeX)/VolumeX;

    Real *AverageSpatialGradient    = (Real*)malloc(sizeof(Real)*ndim*nodeperelem);
    std::fill(AverageSpatialGradient,AverageSpatialGradient+nodeperelem*ndim,0.);
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer i=0; i<ndim; ++i) {
            for (Integer g=0; g<ngauss; ++g) {
                AverageSpatialGradient[a*ndim+i] += SpatialGradient[g*nodeperelem*ndim+a*ndim+i]*detJ[g];
            }
        }
    }
    for (Integer a=0; a<nodeperelem; ++a) {
        for (Integer b=0; b<nodeperelem; ++b) {
            for (Integer i=0; i<ndim; ++i) {
                for (Integer j=0; j<ndim; ++j) {
                    volumetric_stiffness[(a*nvar+i)*nodeperelem*nvar+(b*nvar+j)] += 
                             AverageSpatialGradient[a*ndim+i]*AverageSpatialGradient[b*ndim+j]/VolumeX;
                }
            }
        }
    }
    free(AverageSpatialGradient);
    }
}

#endif //VOLUMETRICSTIFFNESS_H
