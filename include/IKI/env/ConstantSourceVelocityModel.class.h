#ifndef IKI_env_ConstantSourceVelocityModel_H
#define IKI_env_ConstantSourceVelocityModel_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/env/AbstractSourceVelocityModel.interface.h>

namespace IKI { namespace env {
    template <typename T>
    class ConstantSourceVelocityModel final : public AbstractSourceVelocityModel<T> {
    public:
        VectorH<T> at(VectorSp<T> R) const override {
            return v0;
        }

        ConstantSourceVelocityModel(VectorH<T> v0) : v0(v0) { }

    private:
        VectorH<T> v0;
    };
}//env
}//IKI

#endif //IKI_env_ConstantSourceVelocityModel_H