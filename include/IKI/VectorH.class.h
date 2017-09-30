#ifndef IKI_VectorH_CLASS_H
#define IKI_VectorH_CLASS_H

#include <cmath>

#include <IKI/VectorSp.class.h>

namespace IKI {
    template<typename T>
    struct VectorH {
        VectorH(): VectorH(T(0),T(0)) { }
        VectorH(const T& pl_in, const T& pr_in): pl(pl_in), pr(pr_in) { }

        VectorH<T> &operator+=(VectorH<T> const &vec) {
            this->pl += vec.pl;
            this->pr += vec.pr;
            return *this;
        }

        VectorH<T> &operator-=(VectorH<T> const &vec) {
            this->pl -= vec.pl;
            this->pr -= vec.pr;
            return *this;
        }

        VectorH<T> &operator*=(T const &c) {
            this->pl *= c;
            this->pr *= c;
            return *this;
        }

        VectorH<T> &operator/=(T const &c) {
            this->pl /= c;
            this->pr /= c;
            return *this;
        }

        T pl, pr;
    };

    template<typename T>
    VectorH<T> operator+(VectorH<T> const &lha, VectorH<T> const &rha) {
        VectorH<T> tmp(lha); tmp += rha;
        return tmp;
    }

    template<typename T>
    VectorH<T> operator-(VectorH<T> const &lha, VectorH<T> const &rha) {
        VectorH<T> tmp(lha); tmp -= rha;
        return tmp;
    }

    template<typename T>
    VectorH<T> operator*(VectorH<T> const &vec, T const &c) {
        VectorH<T> tmp(vec); tmp *= c;
        return tmp;
    }

    template<typename T>
    VectorH<T> operator*(T const &c,VectorH<T> const &vec) {
        VectorH<T> tmp(vec); tmp *= c;
        return tmp;
    }

    template<typename T>
    VectorH<T> operator/(VectorH<T> const &vec, T const &c) {
        VectorH<T> tmp(vec); tmp /= c;
        return tmp;
    }

    template<typename T>
    VectorH<T> operator/(T const &c,VectorH<T> const &vec) {
        VectorH<T> tmp(vec); tmp /= c;
        return tmp;
    }

    template<typename T>
    T operator*(VectorH<T> const &lha, VectorH<T> const &rha) {
        return lha.pl*rha.pl + lha.pr*rha.pr;
    }

    template<typename T>
    T norm(VectorH<T> const &vec) {
        return sqrt(vec*vec);
    }

    template <typename T>
    VectorH<T> direction_of(VectorH<T> const &of) {
        return of/norm(of);
    }
    
    template <typename T>
    VectorH<T> projection_of_on(VectorSp<T> const &of, VectorSp<T> const &on) {
        T pl = of*on/norm(on);
        T pr = of*of - pl*pl;
        pr = pr > static_cast<T>(0) ? std::sqrt(pr) : static_cast<T>(0);
        return VectorH<T>(pl,pr);
    }
} //IKI

#endif //IKI_VectorH_CLASS_H
