#pragma once

#include "include/BaseAnalysis.h"

using namespace std;

class DT_CMS : public BaseAnalysis {
    public:

    ~DT_CMS();
    void init();

    void Execute(std::mt19937 &engine);
    
    void Finalise() { return; } // maybe print cutflow?
    DT_CMS() { this->setup(); this->analysisname="DT_CMS"; this->init(); };
    
    
};