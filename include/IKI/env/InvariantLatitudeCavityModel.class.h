#ifndef IKI_env_InvariantLatitudeCavityModel_H
#define IKI_env_InvariantLatitudeCavityModel_H

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/env/AbstractCavityModel.interface.h>

namespace IKI { namespace env {
    template <typename T>
    class InvariantLatitudeCavityModel final: public AbstractCavityModel<T> {
    public:
        virtual T at(VectorSp<T> R) const override {
            T L0 = invariant_latitude_calculation(R0);
            T L = invariant_latitude_calculation(R);
            return std::exp(-(L-L0)*(L-L0)/cavity_IL_size_sqr -(R.phi-R0.phi)*(R.phi-R0.phi)/cavity_phi_size_sqr);
        }

        InvariantLatitudeCavityModel(
            VectorSp<T> R0
            ,T cavity_IL_size
            ,T cavity_phi_size
        )
        : R0(R0), cavity_IL_size_sqr(cavity_IL_size*cavity_IL_size), cavity_phi_size_sqr(cavity_phi_size*cavity_phi_size) {
        }

    private:
        T invariant_latitude_calculation(VectorSp<T> R) const {
            return acos(sin(R.th)*sqrt(1./R.r));
        }

        VectorSp<T> R0;
        T cavity_IL_size_sqr, cavity_phi_size_sqr;
    };
}//env
}//IKI

#endif//IKI_env_InvariantLatitudeCavityModel_H
