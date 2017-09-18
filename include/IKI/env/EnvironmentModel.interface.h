#ifndef IKI_env_EnvironmentModel_H
#define IKI_env_EnvironmentModel_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>

namespace IKI { namespace env {
    template<typename T>
    struct EnvironmentModel {
        virtual VectorSp<T> magnetic_field(VectorSp<T> R) const =0;
        
        virtual VectorH<T> source_velocity(VectorSp<T> R) const =0;
        
        virtual T density_cold(VectorSp<T> R) const =0;
        virtual T density_source(VectorSp<T> R) const =0;

        virtual T omega_plasma_cold(VectorSp<T> R) const =0;
        virtual T omega_plasma_source(VectorSp<T> R) const =0;
        virtual T omega_cyclotron_cold(VectorSp<T> R) const =0;

        virtual ~EnvironmentModel() noexcept { };
    };
}//env
}//IKI

#endif //IKI_EnvironmentModel_H
