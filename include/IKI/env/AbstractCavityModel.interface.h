#ifndef IKI_env_AbstractCavityModel_CLASS_H
#define IKI_env_AbstractCavityModel_CLASS_H

#include <IKI/VectorSp.class.h>

namespace IKI { namespace env {
    template <typename T>
    struct AbstractCavityModel {
    public:
        virtual T at(VectorSp<T> R) const =0;
        virtual ~AbstractCavityModel() noexcept { }
    };
}//env
}//IKI

#endif // IKI_env_AbstractCavityModel_CLASS_H
