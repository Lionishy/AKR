#ifndef IKI_go_AbstractDispersionRelation_H
#define IKI_go_AbstractDispersionRelation_H

#include <IKI/VectorH.class.h>
#include <IKI/VectorSp.class.h>
#include <IKI/math/VectorFunction.h>

namespace IKI { namespace go {
    template <typename T>
    struct AbstractDispersionRelation {
        virtual FVector<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const =0;
        virtual FVector<T> at(VectorSp<T> R, VectorH<T> K, FVector<T> w) const =0;

        virtual ~AbstractDispersionRelation() noexcept { }
    };
}//go
}//IKI

#endif //IKI_go_AbstractDispersionRelation_H
