#ifndef IKI_go_AbstractDispersionRelationCorrector_H
#define IKI_go_AbstractDispersionRelationCorrector_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>

namespace IKI { namespace go {
    template <typename T>
    struct AbstractDispersionRelationCorrector {
        virtual int correct(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out) const =0;

        virtual ~AbstractDispersionRelationCorrector() noexcept { }
    };    
}//go
}//IKI

#endif //IKI_go_AbstractDispersionRealationSolver
