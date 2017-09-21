#ifndef IKI_go_AbstactStepLogger_H
#define IKI_go_AbstactStepLogger_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>

namespace IKI { namespace go {
    template <typename T>
    struct AbstractStepLogger {
        virtual void log(int res, VectorSp<T> resR, VectorSp<T> resK, FVector<T> resW, T dt, VectorSp<T> vr, VectorSp<T> vk) const =0;

        virtual ~AbstractStepLogger() noexcept { }
    };
}//go
}//IKI

#endif //IKI_go_AbstractStepLogger_H
