#ifndef IKI_SimpleAsciiStepLogger_H
#define IKI_SimpleAsciiStepLogger_H

#include <iostream>

#include <IKI/go/AbstractStepLogger.interface.h>

namespace IKI {
    template <typename T>
    class SimpleAsciiStepLogger final : public go::AbstractStepLogger<T> {
    public:    
        void log(int res, VectorSp<T> resR, VectorSp<T> resK, FVector<T> resW, T dt, VectorSp<T> vr, VectorSp<T> vk) const override {
            if (0 == res) {

                t += dt;
                S += resW[1]*dt*2.*6380./3.*2.*3.14159265;
                if (skip_number == skip_count++) {
                    skip_count = 0;
                    out
                        << t << " " << dt 
                        << " " << resR.r << " " << resR.th << " " << resR.phi
                        << " " << resK.r << " " << resK.th << " " << resK.phi
                        << " " << vr.r << " " << vr.th << " " << vr.phi
                        << " " << resW[1] << " " << S 
                        << std::endl
                    ;    
                    
                }    
            }
        }

        SimpleAsciiStepLogger<T>& reset(T t, T S, unsigned int skip_count) {
            this->t = t;
            this->S = S;
            this->skip_count = skip_count;
            return *this;
        }

        SimpleAsciiStepLogger(std::ostream &out, T t, T S, unsigned int skip_number) : out(out), t(t), S(S), skip_number(skip_number), skip_count(0U){ }

        SimpleAsciiStepLogger(SimpleAsciiStepLogger<T> const &src) = delete;
        SimpleAsciiStepLogger<T>& operator=(SimpleAsciiStepLogger<T> const &src) = delete;
        SimpleAsciiStepLogger(SimpleAsciiStepLogger<T> &&src) = default;
        SimpleAsciiStepLogger<T>& operator=(SimpleAsciiStepLogger<T> &&src) = default;


    private:
        std::ostream &out;
        mutable T t, S;
        mutable unsigned int skip_count;
        unsigned int const skip_number;
         
    };
}//IKI

#endif //IKI_SimpleAsciiStepLogger_H
