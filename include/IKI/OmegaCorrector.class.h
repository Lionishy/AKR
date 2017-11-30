#ifndef IKI_OmegaCorrector_H
#define IKI_OmegaCorrector_H

#include <memory>

#include <IKI/math/BroydnSolver.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>

namespace IKI {
    template <typename T>
    class OmegaCorrector final : public go::AbstractDispersionRelationCorrector<T> {
    public:
        int correct(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out) const override {
            auto OmegaClosure =
                [&R,&K,this] (FVector<double> x) -> FVector<double> {
                    return dr->at(R,K,x);
                };

            FVector<double> outX, outF;
            int res =
                SolveBroydn(
                      OmegaClosure
                    , w
                    , tolerance 
                    , max_iter 
                    , outX, outF
                    , dX, dY
                );

            if (!res) {
                R_out = R;
                K_out = K;
                w_out = outX;
            }
            return res;
        }

        int correct(VectorSp<T> R, VectorH<T> K, FVector<T> w, VectorSp<T> &R_out, VectorH<T> &K_out, FVector<T> &w_out) const {
            auto OmegaClosure =
                [&R,&K,this] (FVector<double> x) -> FVector<double> {
                    return dr->at(R,K,x);
                };

            FVector<double> outX, outF;
            int res =
                SolveBroydn(
                      OmegaClosure
                    , w
                    , tolerance 
                    , max_iter 
                    , outX, outF
                    , dX, dY
                );

            if (!res) {
                R_out = R;
                K_out = K;
                w_out = outX;
            }
            return res;
        }

        OmegaCorrector(
            std::shared_ptr<go::AbstractDispersionRelation<T> const> dr
            , T tolerance, T dX, T dY, unsigned int max_iter
        ) : dr(dr), tolerance(tolerance), dX(dX), dY(dY), max_iter(max_iter) {
        }    

    private:
        T tolerance, dX, dY;
        unsigned int max_iter;
        std::shared_ptr<go::AbstractDispersionRelation<T> const> dr;
    };
}

#endif //IKI_OmegaCorrector_H

