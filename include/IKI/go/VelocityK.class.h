#ifndef IKI_go_VelocityK_H
#define IKI_go_VelocityK_H

#include <memory>
#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>
#include <IKI/go/AbstractVelocityK.interface.h>

namespace IKI { namespace go {
    template <typename T>
    class VelocityK final : public AbstractVelocityK<T> {
    public:
        VectorSp<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const override {
            auto dr_dr = VectorSp<T>(
                  CentralDifferenceDerivative([&R,&K,&w,this] (T r) -> T { return (dr->at(VectorSp<T>(r,R.th,R.phi),K,w))[0]; }, R.r, dR.r)
                , CentralDifferenceDerivative([&R,&K,&w,this] (T th) -> T { return (dr->at(VectorSp<T>(R.r,th,R.phi),K,w))[0]; }, R.th, dR.th)
                , CentralDifferenceDerivative([&R,&K,&w,this] (T phi) -> T { return (dr->at(VectorSp<T>(R.r,R.th,phi),K,w))[0]; }, R.phi, dR.phi)
            );

            auto dD_dw = CentralDifferenceDerivative([&R,&K,&w,this] (T wr) -> T { return (dr->at(R,K,FVector<T>({wr,w[1]})))[0]; }, w[0], dw);
           
            return VectorSp<T>(
                  dr_dr.r/dD_dw
                , dr_dr.th/dD_dw/R.r
                , dr_dr.phi/dD_dw/R.r/sin(R.th)       
            );        
        }

    public:
        VelocityK(VectorSp<T> dR, T dw, std::shared_ptr<AbstractDispersionRelation<T> const> dr) : dR(dR), dw(dw), dr(dr) {
        }

    private:
        VectorSp<T> dR;
        T dw;
        std::shared_ptr<AbstractDispersionRelation<T> const> dr;
    };
}//go
}//IKI

#endif //IKI_go_VelocityK_H
