#ifndef IKI_KStepper_H
#define IKI_KStepper_H

#include <memory>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractDispersionRelationCorrector.interface.h>

namespace IKI {
    template <typename T>
    class KStepper final {
    public:    
        int find(VectorSp<T> R, VectorSp<T> Kbegin, VectorSp<T> Kend, FVector<T> w, VectorSp<T> &outR, VectorSp<T> &outK, FVector<T> &outW) const {
            auto dK = direction_of(Kend - Kbegin)*step_norm;
            
            while (norm(Kend - Kbegin) > step_norm) {
                int res = corrector->correct(R,Kbegin+dK,w,R,Kbegin,w);
                if ( res )
                    return res;
            }
            
            return corrector->correct(R,Kend,w,outR,outK,outW);
        }

        KStepper(std::shared_ptr<go::AbstractDispersionRelationCorrector<T> const> corrector, T step_norm) : corrector(corrector), step_norm(step_norm) { }
    
    private:    
        std::shared_ptr<go::AbstractDispersionRelationCorrector<T> const> corrector;
        T step_norm;

    };
}//IKI

#endif //IKI_KrStepper_H
