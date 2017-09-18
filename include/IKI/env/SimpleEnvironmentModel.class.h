#ifndef IKI_env_SimpleEnvironmentModel_CLASS_H
#define IKI_env_SimpleEnvironmentModel_CLASS_H

#include <memory>
#include <cmath>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/env/EnvironmentModel.interface.h>
#include <IKI/env/AbstractMagneticFieldModel.interface.h>
#include <IKI/env/AbstractCavityModel.interface.h>
#include <IKI/env/AbstractSourceVelocityModel.interface.h>

namespace IKI { namespace env {
    template <typename T>
    class SimpleEnvironmentModel final: public EnvironmentModel<T> {
    public:    
        VectorSp<T> magnetic_field(VectorSp<T> R) const override {
            return this->field->at(R);
        }
        
        VectorH<T> source_velocity(VectorSp<T> R) const override {
            return this->velocity->at(R); 
        }

        T density_cold(VectorSp<T> R) const override {
            return (R0.r*R0.r)/(R.r*R.r)*(1. - this->cavity->at(R));
        }

        T density_source(VectorSp<T> R) const override {
            return source_density_coeff*(R0.r*R0.r)/(R.r*R.r)*this->cavity->at(R);
        }

        T omega_plasma_cold(VectorSp<T> R) const override {
            return omega_pc0*std::sqrt(this->density_cold(R));
        }

        T omega_plasma_source(VectorSp<T> R) const override {
            return omega_pc0*std::sqrt(this->density_source(R));
        }


        T omega_cyclotron_cold(VectorSp<T> R) const override {
            return omega_cc0*norm(this->field->at(R))/norm(this->field->at(R0));
        }

        SimpleEnvironmentModel(
            T omega_pc0, T omega_cc0, T source_density_coeff
            , VectorSp<T> R0
            , std::shared_ptr<AbstractMagneticFieldModel<T> const> field
            , std::shared_ptr<AbstractCavityModel<T> const> cavity
            , std::shared_ptr<AbstractSourceVelocityModel<T> const> velocity
        ): omega_pc0(omega_pc0), omega_cc0(omega_cc0), source_density_coeff(source_density_coeff), R0(R0), field(field), cavity(cavity), velocity(velocity)  {
        }


    private:
        T omega_pc0, omega_cc0, source_density_coeff;
        VectorSp<T> R0;
        std::shared_ptr<AbstractMagneticFieldModel<T> const> field;
        std::shared_ptr<AbstractCavityModel<T> const> cavity;
        std::shared_ptr<AbstractSourceVelocityModel<T> const> velocity;
    };
}//env
}//IKI

#endif//IKI_SimpleEnvironmentModel_CLASS_H

