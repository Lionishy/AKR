#ifndef IKI_go_AbstractVelocityK_H
#define IKI_go_AbstractVelocityK_H

#include <IKI/math/VectorFunction.h>
#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>

namespace IKI { namespace go {
    template <typename T>
    struct AbstractVelocityK {
        VectorSp<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const =0;
        virtual ~AbstractVelocityK() noexcept { }
    };   
}//go
}//IKI

#endif //IKI_go_AbstractVelocityK_H

