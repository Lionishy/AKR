#ifndef IKI_GammaKperpCorrector_H
#define IKI_GammaKperpCorrector_H

#include <memory>
#include <cmath>

#include <IKI/env/AbstractMagneticFieldModel.interface.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>
#include <IKI/math/BroydnSolver.h>
#include <IKI/math/VectorFunction.h>

namespace IKI {
    template <typename T>
    class GammaKperpCorrector final : public go::AbstractDispersionRelationCorrector<T> {
    public:
        int correct(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out) const override {
            VectorSp<T> magnetic_field = field->at(R);
            VectorH<T> Kh = projection_of_on(K,magnetic_field);

            auto KprGammaClosure =
                [&R,&Kh,&w,this] (FVector<T> x) -> FVector<T> {
                    return dr->at(R,VectorH<T>(Kh.pl,x[0]),FVector<T>({w[0],x[1]}));
                };

            FVector<T> outX, outF;
            FVector<T> x0({Kh.pr,w[1]});
            int res =
                SolveBroydn(
                    KprGammaClosure
                    , x0
                    , tolerance
                    , max_iter
                    , outX
                    , outF
                    , dX, dY
                ); 
            if (0 == res) {
                R_out = R;
                K_out = restore(K,VectorH<T>(Kh.pl,outX[0]),magnetic_field);
                w_out = FVector<T>({w[0],outX[1]});
            }

            return res;
        }

        GammaKperpCorrector(
              std::shared_ptr<go::AbstractDispersionRelation<T> const> dr
            , std::shared_ptr<env::AbstractMagneticFieldModel<T> const> field       
            , T tolerance, T dX, T dY, unsigned int max_iter
        ) : dr(dr), field(field), tolerance(tolerance), dX(dX), dY(dY), max_iter(max_iter) {
        }    

    private:
        VectorSp<T> restore(VectorSp<T> old_K, VectorH<T> new_Projection, VectorSp<T> magnetic_field) const {
            VectorH<T> old_Projection = projection_of_on(old_K,magnetic_field);
            VectorSp<T> He = direction_of(magnetic_field);
            
            T Kphi = old_K.phi*(new_Projection.pr/old_Projection.pr);
            T Kpr_rth = std::sqrt(new_Projection.pr*new_Projection.pr - Kphi*Kphi);
            T coeff = ((VectorSp<T>(-He.th,He.r,0.)*old_K) < 0. ? -1. : 1.);
            T Kr = old_Projection.pl*He.r - Kpr_rth*coeff*He.th;
            T Kth = old_Projection.pl*He.th + Kpr_rth*coeff*He.r;

            return VectorSp<T>(Kr,Kth,Kphi);
        }
        
    private:
        std::shared_ptr<go::AbstractDispersionRelation<T> const> dr;
        std::shared_ptr<env::AbstractMagneticFieldModel<T> const> field;
        T tolerance, dX, dY;
        unsigned int max_iter;
    };
}//IKI

#endif //IKI_GammaKperpCorrector_H
