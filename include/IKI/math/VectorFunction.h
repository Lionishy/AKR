#ifndef IKI_VectorFunction_H
#define IKI_VectorFunction_H

#include <array>
#include <functional>



namespace IKI {
    template<typename X, unsigned int N>
    using VFunction = std::function<std::array<X,N> (std::array<X,N> const &)>;
    
    template<unsigned int P, unsigned int N, typename X>
    std::function<X(std::array<X,N> const &)> MakeSlice(VFunction<X,N> vf) {
        return [vf] (std::array<X,N> const &x) -> X { return vf(x)[P]; }; 
    }
}

#endif //IKI_VectorFunction_H
