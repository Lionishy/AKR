#ifndef IKI_SubsteppingDecorator_H
#define IKI_SubsteppingDecorator_H

#include <memory>

#include <IKI/go/AbstractStep.interface.h>

namespace IKI {
    template <typename T, typename Col>
    class SubsteppingDecorator final : public go::AbstractStep<T> {
    public:
        int step(VectorSp<T> R, VectorSp<T> K, FVector<T> w, VectorSp<T> &R_out, VectorSp<T> &K_out, FVector<T> &w_out, T &dt) const {
            int res;
            for (auto dt_iterator = dt_begin; dt_iterator != dt_end; ++dt_iterator) {
                dt = *dt_iterator;
                if (!(res = inner_step->step(R,K,w,R_out,K_out,w_out,dt)))
                    break;
            }

            return res;
        }

        SubsteppingDecorator(std::shared_ptr<go::AbstractStep<T> const> inner_step, Col collection) : inner_step(inner_step), dt_begin(std::begin(collection)), dt_end(std::end(collection)) {
        }

    private:
        std::shared_ptr<go::AbstractStep<T> const> inner_step;
        Col collection;
        typename Col::iterator const dt_begin, dt_end;
    };
}//IKI

#endif //IKI_SubsteppingDecorator_H
