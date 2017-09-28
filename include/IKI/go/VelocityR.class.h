#ifndef IKI_go_VelocityR_H
#define IKI_go_VelocityR_H

//DEBUG!!!
    #include <iostream>

#include <memory>
#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>
#include <IKI/go/AbstractVelocityR.interface.h>

namespace IKI { namespace go {
    template <typename T>
    class VelocityR final : public AbstractVelocityR<T> {
    public:
        VectorSp<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const override {
            auto dr_dk = VectorSp<T>(
                  CentralDifferenceDerivative([&R,&K,&w,this] (T kr) -> T { return (dr->at(R,VectorSp<T>(kr,K.th,K.phi),w))[0]; }, K.r, dK.r)
                , CentralDifferenceDerivative([&R,&K,&w,this] (T kth) -> T { return (dr->at(R,VectorSp<T>(K.r,kth,K.phi),w))[0]; }, K.th, dK.th)
                , CentralDifferenceDerivative([&R,&K,&w,this] (T kphi) -> T { return (dr->at(R,VectorSp<T>(K.r,K.th,kphi),w))[0]; }, K.phi, dK.phi)
            );

            auto dr_dw = CentralDifferenceDerivative([&R,&K,&w,this] (T wr) -> T { return (dr->at(R,K,FVector<T>({wr,w[1]})))[0]; }, w[0], dw);
            
            return -1.*VectorSp<T>(
                  dr_dk.r/dr_dw
                , dr_dk.th/dr_dw/R.r
                , dr_dk.phi/dr_dw/R.r/sin(R.th)       
            );        
        }

    public:
        VelocityR(VectorSp<T> dK, T dw, std::shared_ptr<AbstractDispersionRelation<T> const> dr) : dK(dK), dw(dw), dr(dr) {
        }

    private:
        VectorSp<T> dK;
        T dw;
        std::shared_ptr<AbstractDispersionRelation<T> const> dr;
    };
}//go
}//IKI


#endif //IKI_go_VelocityR_H
