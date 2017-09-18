#ifndef IKI_env_AbstractMagneticFieldModel_H
#define IKI_env_AbstractMagneticFieldModel_H

#include <IKI/VectorSp.class.h>

namespace IKI { namespace env {
    template <typename T>
    struct AbstractMagneticFieldModel {
        virtual VectorSp<T> at(VectorSp<T> R) const =0;
        virtual ~AbstractMagneticFieldModel() { }
    };
}//env
}//IKI

#endif // IKI_env_AbstractMagneticFieldModel_H
