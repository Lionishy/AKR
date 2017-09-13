#ifndef IKI_Derivative_H
#define IKI_Derivative_H

#include <functional>

namespace IKI {
    template<typename F, typename X>
    std::function<X(X)> CreateForwardDifference(F f, X difference) {
        return [f,difference] (X x) -> X { return f(x+difference) - f(x); };
    }

    template<typename F, typename X>
    std::function<X(X)> CreateBackwardDifference(F f, X difference) {
        return [f,difference] (X x) -> X { return f(x) - f(x-difference); };
    }

    template<typename F, typename X>
    std::function<X(X)> CreateCentralDifference(F f, X difference) {
        return [f,difference] (X x) -> X { return f(x + difference/2) - f(x - difference/2); };
    }

    template<typename DF, typename X>
    std::function<X(X)> CreateDerivative(DF df, X difference) {
        return [df,difference] (X x) -> X { return df(x)/difference; };
    }
}//IKI

#endif
