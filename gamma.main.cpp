#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <array>
#include <vector>
#include <utility>
#include <limits>
#include <string>
#include <sstream>

#include <IKI/VectorSp.class.h>
#include <IKI/VectorH.class.h>
#include <IKI/math/VectorFunction.h>
#include <IKI/OmegaCorrector.class.h>
#include <IKI/env/DipoleFieldModel.class.h>
#include <IKI/env/InvariantLatitudeCavityModel.class.h>
#include <IKI/env/AdiabaticInvariantSourceVelocityModel.class.h>
#include <IKI/env/SimpleEnvironmentModel.class.h>
#include <IKI/dsp/AKRDispersionRelation.class.h>

using namespace IKI;

int main(int argc, char ** argv) {
    //r th phi
    VectorSp<double> R0(1.975,0.44928364,0.);
    //Vparal Vperp
    VectorH<double> V0(0.1,0.12241449);

    auto magnetic_field = std::make_shared<env::DipoleFieldModel<double>>();
    auto cavity = std::make_shared<env::InvariantLatitudeCavityModel<double>>(R0,0.005,0.23); 
    auto source_velocity = std::make_shared<env::AdiabaticInvariantSourceVelocityModel<double>>(R0,V0,magnetic_field);
    auto environment = 
        std::make_shared<env::SimpleEnvironmentModel<double>>(
            0.1,1.0,0.2
            , R0
            , magnetic_field    
            , cavity
            , source_velocity
        );

    auto AKR_dispersion_relation = std::make_shared<dsp::AKRDispersionRelation<double>>(environment);
    
    //start
    VectorSp<double> R(1.97471333,0.44924864,0.00004898);
    VectorSp<double> K(0.0,0.0,0.8);
    FVector<double>  w0({0.989513,0.0050196});

    //reflex
    /*VectorSp<double> R(1.95338729,0.44647380,0.01008005);
    VectorSp<double> K(0.40909074,0.05076166,0.85359524);
    FVector<double> w0({0.989513,0.00536094});*/


    auto omega_corrector = std::make_shared<OmegaCorrector<double>>(AKR_dispersion_relation,1.e-9,1.e-7,1.e-7,1000);
    
    
    double KhStep = 1.e-3;
    //up
    {
        VectorH<double> Kh = projection_of_on(K,magnetic_field->at(R));
        int res = 0;
        FVector<double> w = w0;

        std::ofstream gamma_up_out("../data/gamma.up.txt");
        do {
            res = omega_corrector->correct(R,Kh,w,R,Kh,w);
            if (!res)
                gamma_up_out<< Kh.pl << " " << w[1] << std::endl;
            Kh.pl += KhStep;
        } while (w[1] > 1.e-4 && !res);
    }

    //down
    {
        VectorH<double> Kh = projection_of_on(K,magnetic_field->at(R));
        int res = 0;
        FVector<double> w = w0;

        std::ofstream gamma_dn_out("../data/gamma.dn.txt");
        do {
            res = omega_corrector->correct(R,Kh,w,R,Kh,w);
            if (!res)
                gamma_dn_out<< Kh.pl << " " << w[1] << std::endl;
            Kh.pl -= KhStep;
        } while (w[1] > 1.e-4 && !res);
    }

}