

#include <IKI/math/VectorFunction.h>
#include <IKI/math/BroydnSolver.h>
#include <IKI/env/DipoleFieldModel.class.h>
#include <IKI/env/InvariantLatitudeCavityModel.class.h>
#include <IKI/env/AdiabaticInvariantSourceVelocityModel.class.h>
#include <IKI/env/ConstantSourceVelocityModel.class.h>
#include <IKI/env/SimpleEnvironmentModel.class.h>
#include <IKI/dsp/AKRDispersionRelation.class.h>
#include <IKI/go/VelocityK.class.h>
#include <IKI/go/VelocityR.class.h>
#include <IKI/go/FieldGuidedVelocityR.class.h>
#include <IKI/go/PredictorStep.class.h>
#include <IKI/go/StaggeredStep.class.h>
#include <IKI/GammaKpeprCorrector.class.h>
#include <IKI/OmegaCorrector.class.h>
#include <IKI/SimpleAsciiStepLogger.class.h>
#include <IKI/SubsteppingDecorator.class.h>
#include <IKI/KStepper.class.h>
#include <IKI/Trajectory.class.h>

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

using namespace IKI;

int main(int argc, char ** argv) {
    //r th phi
    VectorSp<double> R0(1.975,0.44928364,0.);
    //Vparal Vperp
    VectorH<double> V0(0.15,0.12247449);
    
    auto magnetic_field = std::make_shared<env::DipoleFieldModel<double>>();
    auto cavity = std::make_shared<env::InvariantLatitudeCavityModel<double>>(R0,0.005,0.23); 
    auto source_velocity = 
       //std::make_shared<env::AdiabaticInvariantSourceVelocityModel<double>>(R0,V0,magnetic_field);
       std::make_shared<env::ConstantSourceVelocityModel<double>>(V0);
    auto environment = 
        std::make_shared<env::SimpleEnvironmentModel<double>>(
            0.1,1.0,0.2
            , R0
            , magnetic_field    
            , cavity
            , source_velocity
        );

    auto AKR_dispersion_relation = std::make_shared<dsp::AKRDispersionRelation<double>>(environment);

    VectorSp<double> K0(0.,0.,0.9);
    //FVector<double> w0({0.976,0.0063}); //V0 = 0.2
    //FVector<double> w0({0.991,0.006}); //V0 = 0.1
    FVector<double> w0({0.985,0.0061}); //V0 = 0.15
    auto omega_corrector = std::make_shared<OmegaCorrector<double>>(AKR_dispersion_relation,1.e-9,1.e-7,1.e-7,1000);
    if (omega_corrector->correct(R0,K0,w0,R0,K0,w0)) {
        std::cout << "Error! Poor initial omega guess..." << std::endl;
        return 0;
    }

    if (4 == argc)
    { 
        std::stringstream k_string; 
        k_string << argv[1] << " " << argv[2] << " " << argv[3];

        VectorSp<double> Kend;
        k_string >> Kend.r >> Kend.th >> Kend.phi;

        std::cout << Kend.r << " " << Kend.th << " " << Kend.phi << std::endl;
        
        KStepper<double> kr_stepper(omega_corrector,0.001);
        {
            int res = kr_stepper.find(R0,K0,Kend,w0,R0,K0,w0);
            if (res) {
                std::cout << "Error: " << res << " " << K0.r;
                return 0;
            }

            std::cout << "R: " << R0.r << " " << R0.th << " " << R0.phi << std::endl;
            std::cout << "K: " << K0.r << " " << K0.th << " " << K0.phi << std::endl;
            std::cout << "Mr: " << magnetic_field->at(R0).r << std::endl; 
            VectorH<double> Kh = projection_of_on(K0,magnetic_field->at(R0));
            std::cout << "Kh: " << Kh.pl << " " << Kh.pr << std::endl;

            std::cout << w0[0] << " " << w0[1] << std::endl;
            std::cout << "Do you want to proceed? " << std::flush;
            {
                //std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
                std::string ans;
                std::cin >> ans;
                if ("n" == ans || "N" == ans || "no" == ans || "No" == ans) 
                   return 0;
            }       
        }
    }


    std::cout << w0[0] << " " << w0[1] << std::endl;

    //auto default_velocity_R = std::make_shared<go::VelocityR<double>>(VectorSp<double>(1.e-4,1.e-4,1.e-4),1.e-5,AKR_dispersion_relation);
    
    auto field_guided_velocity_R = std::make_shared<go::FieldGuidedVelocityR<double>>(VectorH<double>(1.e-5,1.e-5),1.e-5,AKR_dispersion_relation,magnetic_field);
    auto velocity_K = std::make_shared<go::VelocityK<double>>(VectorSp<double>(1.e-5,1.e-5,1.e-5),1.e-5,AKR_dispersion_relation);
    
    
    auto corrector = std::make_shared<GammaKperpCorrector<double>>(AKR_dispersion_relation,magnetic_field,1.e-8,1.e-4,1.e-4,1000);
    
    

    std::ofstream gain_out("../data/gain.txt");
    gain_out << std::setprecision(8) << std::fixed;
    
    Trajectory<double> trajectory_calculator(field_guided_velocity_R,velocity_K,corrector,2.e-4);
    KStepper<double> kth_stepper(omega_corrector,0.0001);

    bool proceed = true;
    while (proceed) {
        auto logger = std::make_shared<go::EmptyStepLogger<double>>();
        
        double S = 0, dt = 0.1e-4, t = 0;
        VectorSp<double> R, K;
        FVector<double> w;
        int res = trajectory_calculator.trajectory(R0,K0,w0,logger,dt,t,R,K,w,S);

        gain_out << K0.r << " " << K0.th << " " << K0.phi << " " << S << " " << (res?0:1) << std::endl;
        std::cout << K0.th << std::endl;
        {
            VectorSp<double> Kend(K0.r,K0.th-0.01,K0.phi);
            int res = kth_stepper.find(R0,K0,Kend,w0,R0,K0,w0);
            if (w0[1] < 1.e-4 || res) {
                std::cout << "------------------ END ---------------------" << std::endl;
                break;
            }
        }
    }
    return 0;
}
