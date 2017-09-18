#ifndef IKI_env_AdiabaticInvariantSourceVelocityModel_H
#define IKI_env_AdiabaticInvariantSourceVelocityModel_H

#include <memory>
#include <cmath>
#include <cassert>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/env/AbstractMagneticFieldModel.interface.h>
#include <IKI/env/AbstractSourceVelocityModel.interface.h>

namespace IKI { namespace env {
    template <typename T>
    class AdiabaticInvariantSourceVelocityModel final : public  AbstractSourceVelocityModel<T> {
    public:
        VectorH<T> at(VectorSp<T> R) const override {
            VectorSp<T> H = magnetic_field->at(R);
            T v_perp_sqr = v0.pr*v0.pr*norm(magnetic_field->at(R))/norm(magnetic_field->at(R0));
            T v_paral_sqr = v0.pr*v0.pr + v0.pl*v0.pl - v_perp_sqr;

            assert(v_paral_sqr > 0);
            return VectorH<T>(sqrt(v_paral_sqr),sqrt(v_perp_sqr));
        }
       
        AdiabaticInvariantSourceVelocityModel(VectorSp<T> R0, VectorH<T> v0, std::shared_ptr<AbstractMagneticFieldModel<T> const> magnetic_field) : v0(v0), magnetic_field(magnetic_field) {
        }
       
    private:
        VectorSp<T> R0;
        VectorH<T> v0;
        std::shared_ptr<AbstractMagneticFieldModel<T> const> magnetic_field;
    };
}//env
}//IKI

#endif //IKI_env_AdiabaticInvariantSourceVelocityModel_H
