#ifndef IKI_go_Step_H
#define IKI_go_Step_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>

namespace IKI { namespace go {
    template <typename T>
    struct AbstractStep {
        virtual int step(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out, T &dt) const =0;
        
        virtual ~AbstractStep() noexcept { }
    };
}//go
}//IKI

#endif //IKI_go_Step_H
