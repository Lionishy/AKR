#ifndef IKI_VectorFunction_H
#define IKI_VectorFunction_H

#include <array>
#include <functional>
#include <algorithm>

namespace IKI {
    
    /** Vector **/
    template <typename X>
    using FVector = std::array<X,2>;    

    /** template<typename X>
    using VFunction = std::function<std::array<X,N> (std::array<X,N> const &)>;

    template<typename X, unsigned int N>
    using SFunction = std::function<X(std::array<X,N> const &)>;
    
    template<unsigned int P, unsigned int N, typename X>
    SFunction MakeSlice(VFunction<X,N> vf) {
        return [vf] (std::array<X,N> const &x) -> X { return vf(x)[P]; }; 
    }**/

    template <typename X>
    FVector<X> operator+(FVector<X> const &lha, FVector<X> const &rha) {
        return std::array<X,2>({lha[0]+rha[0],lha[1]+rha[1]});
    }

    template <typename X>
    FVector<X> operator-(FVector<X> const &lha, FVector<X> const &rha) {
        return std::array<X,2>({lha[0]-rha[0],lha[1]-rha[1]});
    }

    template <typename X>
    FVector<X> operator*(FVector<X> const &lha, X x) {
        return std::array<X,2>({lha[0]*x,lha[1]*x});
    }
    
    template <typename X>
    FVector<X> operator*(X x,FVector<X> const &lha) {
        return std::array<X,2>({x*lha[0],x*lha[1]});
    }

    template <typename X>
    X operator*(FVector<X> const &lha, FVector<X> const &rha) {
        return lha[0]*rha[0] + rha[1]*lha[1];
    }

    template <typename X>
    FVector<X> operator/(FVector<X> const &lha, X x) {
        return std::array<X,2>({lha[0]/x,lha[1]/x});
    }

    template <typename X>
    X norm(FVector<X> const &vec) {
        return vec*vec;
    }

    template <typename F, typename X>
    X CentralDifferenceDerivative(F f, X x0, X diff) {
        return (f(x0 + diff/static_cast<X>(2)) - f(x0 - diff/static_cast<X>(2)))/diff;
    }

    template <typename F, typename X>
    FVector<X> CentralDifferenceDerivative(F f, FVector<X> x0, FVector<X> diff) {
        return (f(x0+diff/static_cast<X>(2)) - f(x0-diff/static_cast<X>(2)))/norm(diff);
    }

    /** Matrix **/
    template <typename X>
    using FMatrix = std::array<X,4>;

    template <typename X>
    FMatrix<X> operator+(FMatrix<X> const &lha, FMatrix<X> const &rha) {
        return std::array<X,4>({lha[0]+rha[0],lha[1]+rha[1],lha[2]+rha[2],lha[3]+rha[3]});
    }

    template <typename X>
    FMatrix<X> operator-(FMatrix<X> const &lha, FMatrix<X> const &rha) {
        return std::array<X,4>({lha[0]-rha[0],lha[1]-rha[1],lha[2]-rha[2],lha[3]-rha[3]});
    }

    template <typename X>
    FMatrix<X> operator*(FMatrix<X> const &lha, FMatrix<X> const &rha) {
        return std::array<X,4>({
            lha[0]*rha[0] + lha[1]*rha[2]
            ,lha[0]*rha[1] + lha[1]*rha[3]
            ,lha[2]*rha[0] + lha[3]*rha[2]
            ,lha[2]*rha[1] + lha[3]*rha[3]
        });
    }

    template <typename X>
    FVector<X> make_vector(FVector<X> const &lha, FMatrix<X> const &rha) {
        return std::array<X,2>({lha[0]*rha[0] + lha[1]*rha[2], lha[0]*rha[1] + lha[1]*rha[3]});
    }

    template <typename X>
    FVector<X> make_vector(FMatrix<X> const &rha, FVector<X> const &lha) {
        return std::array<X,2>({rha[0]*lha[0] + rha[1]*lha[1], rha[2]*lha[0] + rha[3]*lha[1]});
    }

    template <typename X>
    FMatrix<X> operator*(FMatrix<X> const &m, X c) {
        return FMatrix<X>({m[0]*c,m[1]*c,m[2]*c,m[3]*c});
    }

    template <typename X>
    FMatrix<X> operator*(X c,FMatrix<X> const &m) {
        return std::array<X,2>({c*m[0],c*m[1],c*m[2],c*m[3]});
    }

    template <typename X>
    FMatrix<X> make_matrix(FVector<X> const &lha, FVector<X> const &rha) {
        return std::array<X,4>({
              lha[0]*rha[0]
            , lha[0]*rha[1]
            , lha[1]*rha[0]
            , lha[1]*rha[1]
        });
    }
    
    template <typename X>
    X det(FMatrix<X> const &m) {
        return m[0]*m[3] - m[1]*m[2];
    }    
}//IKI

#endif //IKI_VectorFunction_H
