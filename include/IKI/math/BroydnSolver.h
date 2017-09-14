#ifndef IKI_BroydnSolver_H
#define IKI_BroydnSolver_H

#include <array>
#include <cmath>
#include <IKI/math/Derivative.h>
#include <IKI/math/VectorFunction.h>

namespace IKI {
    
    template <typename X>
    using Jacobian = FMatrix<X>;

    template <typename F, typename X>
    Jacobian<X> CalculateJacobian(F f, FVector<X> x0, X dX, X dY) {
        auto df_dx = CentralDifferenceDerivative(f,x0,FVector<X>({dX,0}));
        auto df_dy = CentralDifferenceDerivative(f,x0,FVector<X>({0,dY}));
        return std::array<X,4>({df_dx[0],df_dy[0],df_dx[1],df_dy[1]});
    }

    template <typename X>
    FVector<X> SolveJacobian(Jacobian<X> J, X detOver, FVector<X> x, FVector<X> b) {
        auto det1 = det(FMatrix<X>({b[0],J[1],b[1],J[3]}));
        auto det2 = det(FMatrix<X>({J[0],b[0],J[2],b[1]}));
        
        return FVector<X>({det1*detOver+x[0] , det2*detOver+x[1]});
    }

    template <typename X>
    Jacobian<X> CalculateNextJacobian(Jacobian<X> jacob, FVector<X> df, FVector<X> dx) {
        return jacob + make_matrix(df-make_vector(jacob,dx),dx)*(static_cast<X>(1)/norm(dx));
    } 

    enum class SolverStatusBroydn {OK=0, NO_PROGRESS=2, MAX_ITERATIONS=4,SINGULAR_JACOBIAN=8};

    template <typename F, typename X>
    SolverStatusBroydn SolveBroydn(
        F f
        , FVector<X> x0
        , X tolerance
        , unsigned int max_iter 
        , FVector<X> &outX
        , FVector<X> &outF
        , X dX, X dY
    ) {
        
        FVector<X> xPrev = x0, fPrev = f(x0), xNext, fNext;
        auto jacob = CalculateJacobian(f,xPrev,dX,dY);

        auto res = SolverStatusBroydn::OK;
        unsigned int iter_made = 0;
        int jacob_sing = 0, no_f_progress = 0, no_x_progress = 0;
        int progress_failed = 0;

        while (norm(fPrev) > tolerance) {
            if (iter_made == max_iter) {
                res = SolverStatusBroydn::MAX_ITERATIONS;
                break;        
            }
            
            X detOver = 1/det(jacob);
            if (std::isnan(detOver) && std::isinf(detOver)) {
                if (1 == jacob_sing) {
                    res = SolverStatusBroydn::SINGULAR_JACOBIAN;
                    break;
                }
                jacob = CalculateJacobian(f,xPrev,dX,dY);
                ++jacob_sing;
                continue;
            }
            else
                jacob_sing = jacob_sing ? --jacob_sing : jacob_sing;

            xNext = SolveJacobian(jacob,detOver,xPrev,static_cast<X>(-1)*fPrev);
            fNext = f(xNext);
            if (!(norm(fNext)<norm(fPrev))) {
                if (1 == no_f_progress) {
                    res = SolverStatusBroydn::NO_PROGRESS;
                    break;
                }
                jacob = CalculateJacobian(f,xPrev,dX,dY);
                ++no_f_progress;
                continue;
            }
            else
                no_f_progress = no_f_progress ? --no_f_progress : no_f_progress;


            auto dx = xNext - xPrev;
            auto dxOver = 1/(dx*dx);
            if (std::isnan(dxOver) || std::isinf(dxOver)) {
                if (1 == no_x_progress) {
                    res = SolverStatusBroydn::NO_PROGRESS;
                    break;
                }
                jacob = CalculateJacobian(f,xPrev,dX,dY);
                ++no_x_progress;
            }
            else
                no_x_progress = no_x_progress ? --no_x_progress : no_x_progress;

            jacob = CalculateNextJacobian(jacob,fNext-fPrev,dx);
            xPrev = xNext;
            fPrev = fNext;
            ++iter_made;    
        }

        outX = xPrev;
        outF = fPrev;
        return res;
    }

}//IKI
#endif //IKI_BroydnSolver_H 


