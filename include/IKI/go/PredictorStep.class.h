#ifndef IKI_go_PredictorStep_H
#define IKI_go_PredictorStep_H

#include <memory>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractStep.interface.h>
#include <IKI/go/VelocityR.class.h>
#include <IKI/go/VelocityK.class.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>

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
            
            return corrector->correct(new_R,new_K,w,R_out,K_out,w_out);            
        }

        PredictorStep(std::shared_ptr<VelocityR<T> const> VR,std::shared_ptr<VelocityK<T> const> VK,std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector) : VR(VR), VK(VK), corrector(corrector) {
        }

    private:
        std::shared_ptr<VelocityR<T> const> VR;
        std::shared_ptr<VelocityK<T> const> VK;
        std::shared_ptr<AbstractDispersionRelationCorrector<T> const> corrector;
    };
}//go
}//IKI

#endif //IKI_go_PredictorStep_H
