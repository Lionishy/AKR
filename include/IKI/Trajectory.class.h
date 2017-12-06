#ifndef IKI_Trajectory_CLASS_H
#define IKI_Trajectory_CLASS_H

#include <IKI/VectorSp.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/go/AbstractStepLogger.interface.h>
#include <IKI/go/StaggeredStep.class.h>

#include <memory>

namespace IKI {
    template <typename T>
    class Trajectory {
        public:
            int trajectory(VectorSp<T> R, VectorSp<T> K, FVector<T> w, std::shared_ptr<go::AbstractStepLogger<T>> logger, T dt, T &t_out, VectorSp<T>  &R_out, VectorSp<T> &K_out, FVector<T> &w_out, T &S_out) {
                auto staggered_step =
                    go::StaggeredStep<double>(velocity_R, velocity_K, corrector, logger,dt,K);

                int res;
                double t = 0, S = 0;
                while ( w[1] > max_gamma) {
                    if (0 < (res = staggered_step.step(R,K,w,R_out,K_out,w_out,dt)) ) {
                        std::cout << t << " " << dt << " Error: " << res << std::endl;
                        break; 
                    }
                    t += dt;
                    R = R_out;
                    K = K_out;
                    w = w_out;
                    S += w[1]*dt*2.*6380./3.*2.*3.14159265;
                }
                S_out = S;
                return res;
            }

            Trajectory(std::shared_ptr<go::AbstractVelocityR<T>> velocity_R, std::shared_ptr<go::AbstractVelocityK<T>> velocity_K, std::shared_ptr<go::AbstractDispersionRelationCorrector<T>> corrector, T max_gamma): velocity_R(velocity_R), velocity_K(velocity_K), corrector(corrector), max_gamma(max_gamma) {
            }
        private:
            T const max_gamma;
            std::shared_ptr<go::AbstractVelocityR<T>> velocity_R;
            std::shared_ptr<go::AbstractVelocityK<T>> velocity_K;
            std::shared_ptr<go::AbstractDispersionRelationCorrector<T>> corrector;
    };
}//IKI

#endif //IKI_Trajectory_CLASS_H