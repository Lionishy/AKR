#ifndef IKI_env_AbstractSourceVelocityModel_H
#define IKI_env_AbstractSourceVelocityModel_H

#include <IKI/VectorH.class.h>

namespace IKI { namespace env {
    template <typename T>
    struct AbstractSourceVelocityModel {
        virtual VectorH<T> at(VectorSp<T> R) const =0;

        virtual ~AbstractSourceVelocityModel() noexcept { }
    };
}//env
}//IKI

#endif // IKI_env_AbstractSourceVelocityModel_H
