#pragma once

#ifndef HSCP_ATLAS_h
#define HSCP_ATLAS_h

#include "include/BaseAnalysis.h"
#include "include/HEPData.h"


using namespace std;

class HSCP_ATLAS : public BaseAnalysis {
    public:

    ~HSCP_ATLAS();
    void init();

    void Execute(std::mt19937 &engine);

    void Finalise() { return; } 
    HSCP_ATLAS() { this->setup(); this->analysisname="HSCP_ATLAS"; this->init(); };
    

   Efficiency1D EtmissTurnOn;   
   Efficiency2D SingleMuTurnOn;
   Efficiency2D LooseEff;   
   Efficiency2D TightPromotionEff;  

   Efficiency1D MToFFullDet;
   Efficiency1D MToFFullDetErr;
    //std::vector<double> massToF2Cand_min;//  = {150.,350.,575.,800.};
    //std::vector<double> massToF1Cand_min;//  = {175.,375.,600.,825.};
};

#endif
