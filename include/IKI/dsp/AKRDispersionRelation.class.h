#ifndef IKI_dsp_AKRDispersionRelation_H
#define IKI_dsp_AKRDispersionRelation_H

#include <memory>
#include <complex>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/env/EnvironmentModel.interface.h>
#include <IKI/go/AbstractDispersionRelation.interface.h>

namespace IKI { namespace dsp {
    template<typename T>
    class AKRDispersionRelation final : public go::AbstractDispersionRelation<T> {
    public:
        FVector<T> at(VectorSp<T> R, VectorSp<T> K, FVector<T> w) const override {
            return at(R,projection_of_on(K,env_model->magnetic_field(R)),w);       
        }

        FVector<T> at(VectorSp<T> R, VectorH<T> K, FVector<T> w) const override {
            std::complex<T> res = at(R,K,std::complex<T>(w[0],w[1]));
            return FVector<T>({res.real(),res.imag()});
        }

        AKRDispersionRelation(std::shared_ptr<env::EnvironmentModel<T> const> env_model): env_model(env_model) {
        }
    
    private:
        std::complex<T> at(const VectorSp<T> &R, const VectorH<T> &k, const std::complex<T> &w) const {
            std::complex<T> w_sqr = w*w;

            T omega_pc = env_model->omega_plasma_cold(R), omega_cc = env_model->omega_cyclotron_cold(R);
            T omega_pcs = omega_pc*omega_pc, omega_ccs = omega_cc*omega_cc;

            T omega_ps = env_model->omega_plasma_source(R);
            T omega_pss = omega_ps*omega_ps;

            VectorH<T> v  = env_model->source_velocity(R);
            T
                sf    = 1. - 0.5*(v.pl*v.pl+v.pr*v.pr),
                kv_pl = k.pl*v.pl;

            std::complex<T>
                n_pl = k.pl/w,
                n_pr = k.pr/w;

            std::complex<T>
                eps1_c = 1. - (omega_pcs)/(w_sqr-omega_ccs),
                eps2_c = omega_cc/w*omega_pcs/(w_sqr-omega_ccs);

            std::complex<T>
                f1 = w-omega_cc*sf-kv_pl,
                f2 = w+omega_cc*sf-kv_pl;

            std::complex<T>
                F1 = (w-kv_pl)/f1 + 0.5*v.pr*v.pr*(k.pl*k.pl-w*omega_cc)/f1/f1,
                F2 = (w-kv_pl)/f2 + 0.5*v.pr*v.pr*(k.pl*k.pl-w*omega_cc)/f2/f2;

            std::complex<T>
                eps1_s = -0.5*omega_pss/w_sqr*(F1+F2),
                eps2_s =  0.5*omega_pss/w_sqr*(F1-F2);

            std::complex<T>
                eps1 = eps1_c+eps1_s,
                eps2 = eps2_c+eps2_s;

            return ((eps1-n_pl*n_pl)/eps1 - n_pr*n_pr)*((eps1*eps1-eps2*eps2)/eps1 - n_pl*n_pl - n_pr*n_pr) - n_pl*n_pl*eps2/eps1*eps2/eps1;
        }

    private:
        std::shared_ptr<env::EnvironmentModel<T> const> env_model;
    };
}//AKR
}//IKI

#endif //IKI_AKR_DispersionRelation_TEMPLATE
