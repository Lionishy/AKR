#ifndef IKI_env_DipoleFieldModel_H
#define IKI_env_DipoleFieldModel_H

#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/env/AbstractMagneticFieldModel.interface.h>

namespace IKI { namespace env {
    template <typename T>
    class DipoleFieldModel final: public AbstractMagneticFieldModel<T> {
    public:
        VectorSp<T> at(VectorSp<T> R) const override {
            return VectorSp<T>(-std::cos(R.th)/(R.r*R.r*R.r),-0.5*std::sin(R.th)/(R.r*R.r*R.r),0.);
        }

        DipoleFieldModel() {
        }
    };
}//env
}//IKI

#endif //IKI_env_DipoleFieldModel_H
