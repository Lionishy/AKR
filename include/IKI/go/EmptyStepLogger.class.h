#ifndef IKI_go_EmptyStepLogger_H
#define IKI_go_EmptyStepLogger_H

#include <IKI/go/AbstractStepLogger.interface.h>

namespace IKI { namespace go {
    template <typename T>
    class EmptyStepLogger final : public AbstractStepLogger<T> {
    public:
        void log(int res, VectorSp<T> resR, VectorSp<T> resK, FVector<T> resW, T dt, VectorSp<T> vr, VectorSp<T> vk) const override { }
    };
}//go
}//IKI
#endif //IKI_go_EmptyStepLogger_H
