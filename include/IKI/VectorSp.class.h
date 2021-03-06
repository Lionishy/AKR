#ifndef IKI_VectorSp_CLASS_H
#define IKI_VectorSp_CLASS_H

#include <cmath>

namespace IKI {
    template<typename T>
    struct VectorSp {
        VectorSp(): VectorSp(T(0),T(0),T(0)) { }
        VectorSp(const T& r_in, const T& th_in, const T& phi_in): r(r_in), th(th_in), phi(phi_in) { }

        const VectorSp<T> &operator+=(const VectorSp<T> &vec) {
            this->r += vec.r;
            this->th += vec.th;
            this->phi += vec.phi;

            return *this;
        }

        const VectorSp<T> &operator-=(const VectorSp<T> &vec) {
            this->r -= vec.r;
            this->th -= vec.th;
            this->phi -= vec.phi;

            return *this;
        }

        const VectorSp<T> &operator*=(const T &c) {
            this->r *= c;
            this->th *= c;
            this->phi *= c;

            return *this;
        }

        const VectorSp<T> &operator/=(const T &c) {
            this->r /= c;
            this->th /= c;
            this->phi /= c;

            return *this;
        }

        

        T r, th, phi;
    };

    template <typename T>
    VectorSp<T> operator+(const VectorSp<T> &lha, const VectorSp<T> &rha) {
        VectorSp<T> tmp(lha); tmp += rha;
        return tmp;
    }

    template <typename T>
    VectorSp<T> operator-(const VectorSp<T> &lha, const VectorSp<T> &rha) {
        VectorSp<T> tmp(lha); tmp -= rha;
        return tmp;
    }

    template <typename T>
    VectorSp<T> operator*(const T &c, const VectorSp<T> &rha) {
        VectorSp<T> tmp(rha); tmp *= c;
        return tmp;
    }

    template <typename T>
    VectorSp<T> operator*(const VectorSp<T> &lha, const T &c) {
        VectorSp<T> tmp(lha); tmp *= c;
        return tmp;
    }

    template <typename T>
    VectorSp<T> operator/(const VectorSp<T> &lha, const T &c) {
        VectorSp<T> tmp(lha); tmp /= c;
        return tmp;
    }

    template <typename T>
    T operator*(const VectorSp<T> &lha, const VectorSp<T> &rha) {
        return lha.r*rha.r + lha.th*rha.th + lha.phi*rha.phi;
    }

    template <typename T>
    T norm(const VectorSp<T> &vec) {
        return sqrt(vec*vec);
    }
    
    template <typename T>
    VectorSp<T> rotate(VectorSp<T> const &t, VectorSp<T> const &R, VectorSp<T> const &dR) {
        return VectorSp<T>(
              t.r + t.th*dR.th + t.phi*dR.phi*sin(R.th)
            , t.th - t.r*dR.th + t.phi*dR.phi*cos(R.th)
            , t.phi - dR.phi*(t.r*sin(R.th) + t.th*cos(R.th))
        );    
    }

    template <typename T>
    VectorSp<T> direction_of(VectorSp<T> const &of) {
        return of/norm(of);
    }
} //IKI

#endif //IKI_VectorSp_CLASS_H
