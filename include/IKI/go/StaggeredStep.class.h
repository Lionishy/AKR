#ifndef IKI_go_StaggeredStep_H
#define IKI_go_StaggeredStep_H

#include <memory>
#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractStep.interface.h>
#include <IKI/go/AbstractVelocityK.interface.h>
#include <IKI/go/AbstractVelocityR.interface.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>
#include <IKI/go/AbstractStepLogger.interface.h>

namespace IKI { namespace go {
    template <typename T>
    class StaggeredStep final : public AbstractStep<T> {
    public:
        int step(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out, T &dt) const override {
            K = (K-K_prev)*(dt/dt_prev) + K_prev;
            
            auto VR = vr->at(R,K,w);
            auto R_new = R + dt*VR;

            auto VK = vk->at(R_new,K,w);
            auto K_new = K + dt*VK;

            K_new = rotate(K_new,R_new,dt*VR);
            
            int res = corrector->correct(R_new,K_new,w,R_out,K_out,w_out);
            
            
            if (!res) {
                logger->log(res,R_out,K_out,w_out,dt,VR,VK);
                K_prev = K;
                dt_prev = dt;
            }

            return res;
        }

        StaggeredStep(
              std::shared_ptr<AbstractVelocityR<T> const> vr
            , std::shared_ptr<AbstractVelocityK<T> const> vk
            , std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector
            , std::shared_ptr<AbstractStepLogger<T> const> logger
            , T dt_prev
            , VectorSp<T> K_prev
        ) : vr(vr), vk(vk), corrector(corrector), logger(logger), dt_prev(dt_prev), K_prev(K_prev) { }
        
    private:
        std::shared_ptr<AbstractVelocityK<T> const> vk;
        std::shared_ptr<AbstractVelocityR<T> const> vr;
        std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector;
        std::shared_ptr<AbstractStepLogger<T> const> logger;
        mutable T dt_prev;
        mutable VectorSp<T> K_prev;
    };
}//go
}//IKI

#endif //IKI_go_StaggeredStep_H