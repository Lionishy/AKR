#ifndef IKI_go_FieldGuidedVelocityR_H
#define IKI_go_FieldGuidedVelocityR_H

#include <memory>
#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>
#include <IKI/env/AbstractMagneticField.interface.h>

namespace IKI { namespace go {
    template <typename T>
    class FieldGuidedVelocityR final {
    public:
        VectorH<T> at(VectorSp<T> R, VectorH<T> K, FVector<T> w) const {
            T dD_dKpl = CentralDifferenceDerivative([&R,&K,&w,this] (T x) -> T { return dr->at(R,VectorH<T>(x,K.pr),w)[0]; },K.pl,dK.pl);
            T dD_dKpr = CentralDifferenceDerivative([&R,&K,&w,this] (T x) -> T { return dr->at(R,VectorH<T>(K.pl,x),w)[0]; },K.pr,dK.pr);

            T dD_dw = CentralDifferenceDerivative([&R,&K,&w,this] (T x) -> T { return dr->at(R,K,FVector<T>({x,w[1]})); },w[0],dw);

            return static_cast<T>(-1)*VectorH<T>(dD_dKpl/dD_dw,dD_dKpr/dD_dw);
        }

        VectorSp<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const {
            auto H = field->at(R);
            auto guided_velocity = this->at(R,projection_of_on(K,H),w);
            auto Kparallel = (K*H)/(H*H)*H;
            auto Kperp = K - Kparallel;
            auto Velocity = guided_velocity.pl*direction_of(Kparallel) + guided_velocity.pr*direction_of(Kperp);

            Velocity.th /= R.r;
            Velocity.phi /= (R.r*std::sin(R.th));

            return Velocity;
        }


    private:
        std::shared_ptr<AbstractDispersionRelation<T> const> dr;
        std::shared_ptr<env::AbstractMagneticField<T> const> field;
        VectorH<T> dK;
        T dw;
    };
}//go
}//IKI

#endif
