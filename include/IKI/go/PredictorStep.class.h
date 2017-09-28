#ifndef IKI_go_PredictorStep_H
#define IKI_go_PredictorStep_H

#include <memory>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractStep.interface.h>
#include <IKI/go/AbstractVelocityR.interface.h>
#include <IKI/go/AbstractVelocityK.interface.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>
#include <IKI/go/AbstractStepLogger.interface.h>
#include <IKI/go/EmptyStepLogger.class.h>

namespace IKI { namespace go {
    template <typename T>
    class PredictorStep final : public AbstractStep<T> {
    public:
        int step(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out, T &dt) const override {
            auto vr = VR->at(R,K,w);            
            auto vk = VK->at(R,K,w);
            
            auto new_R = R + 0.5*dt*vr;
            auto new_K = K + 0.5*dt*vk;
            new_K = rotate(new_K,new_R,0.5*dt*vr);
            
            vr = VR->at(new_R,new_K,w);
            vk = VK->at(new_R,new_K,w);
            
            new_R = R + dt*vr;
            new_K = K + dt*vk;
            new_K = rotate(new_K,new_R,dt*vr);
            
            int res = corrector->correct(new_R,new_K,w,R_out,K_out,w_out);            
            if (res)
                logger->log(res,new_R,new_K,w,dt,vr,vk);
            else
                logger->log(res,R_out,K_out,w_out,dt,vr,vk);

            return res;
        }

        PredictorStep(std::shared_ptr<AbstractVelocityR<T> const> VR,std::shared_ptr<AbstractVelocityK<T> const> VK,std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector): PredictorStep(VR,VK,corrector,std::make_shared<EmptyStepLogger<T>>()) { }

        PredictorStep(std::shared_ptr<AbstractVelocityR<T> const> VR,std::shared_ptr<AbstractVelocityK<T> const> VK,std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector, std::shared_ptr<AbstractStepLogger<T> const> logger) : VR(VR), VK(VK), corrector(corrector), logger(logger) {
        }
        


    private:
        std::shared_ptr<AbstractVelocityR<T> const> VR;
        std::shared_ptr<AbstractVelocityK<T> const> VK;
        std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector;
        std::shared_ptr<AbstractStepLogger<T> const> logger;
    };
}//go
}//IKI

#endif //IKI_go_PredictorStep_H
