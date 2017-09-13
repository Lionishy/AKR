#ifndef IKI_VectorH_CLASS_H
#define IKI_VectorH_CLASS_H

#include <cmath>

namespace AKR {
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
} //IKI

#endif //IKI_VectorH_CLASS_H
