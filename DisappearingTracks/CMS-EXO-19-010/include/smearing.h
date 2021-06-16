

/* 

File for implementing smearing functions. 

For the analyses in version 1 we don't actually use them, hence most are commented out.

Essentially smearing can be implemented in the same way as in GAMBIT; the functions there can be converted for use here with
little work. 

*/




#pragma once

#include <cassert>
#include <cstdlib>
#include <cstddef>
#include <type_traits>

#include <random>
#include <cfloat> // for DBL_MAX


#include "heputil.h"


using namespace std;






template<typename T1> void filterEfficiency(std::vector<T1*> &vec_t1, double eff_fn(T1* t1), std::mt19937 &engine)
{
  //std::vector<bool> mask(vec_t1.size(),false);
  std::uniform_real_distribution<double> rd(0.0,1.0);
  auto it = vec_t1.begin();
  while(it != vec_t1.end())
  {
    double eff=eff_fn(*it);
    if(rd(engine) > eff)
        {
	        it = vec_t1.erase(it);
        }
    else
      {
	it++;
      }
    
  }


}




/*
inline void smearbtag(bool *btag, double pT, std::mt19937 &engine)
{
  if(!(*btag)) return;

  const std::vector<double> binedges_pt = {0,20.,30.,40.,60.,80.,200.,250.,DBL_MAX};
  const std::vector<double> btageff = {0.0,0.54,0.64,0.7,0.74,0.77,0.75,0.7};
  static HEP::BinnedFn1D<double> _btageff1D(binedges_pt,btageff);

  double teff=_btageff1D.get_at(pT);
  static std::uniform_real_distribution<double> rd(0.0,1.0);


  if(rd(engine) > teff)
  {
    *btag=false;
  }
  
}






double JetEff(HEP::Jet* jet) {

  // "medium" working point of JVT figure 17(b) 1510.03823
  // 
  const std::vector<double> binedges_pt = {0,20.,25.,30.,35.,40.,45.,50.,55.,120.,DBL_MAX};
  const std::vector<double> JVTeff = {0.0,0.86,0.9,0.92,0.93,0.94,0.95,0.95,0.96,1.0}; // jets with pT > 120 don't need to satify JVT
  static HEP::BinnedFn1D<double> _JVTeff1D(binedges_pt,JVTeff);


  if(jet->abseta() > 2.8) return 1.0; // vertex tagger to suppress pileup only used for central jets


  double teff=_JVTeff1D.get_at(jet->pT());
  return teff;

}


double ElectronEff(HEP::Particle* electron) {
        static HEP::BinnedFn2D<double> _elEff2d({{0, 1.5, 2.5, DBL_MAX}}, //< |eta|
                                                     {{0, 10., DBL_MAX}}, //< pT
                                                     {{0., 0.95, // |eta| 0.1-1.5
                                                       0., 0.85, // |eta| 1.5-2.5
                                                       0., 0.}}); // |eta| > 2.5
        double teff=_elEff2d.get_at(electron->abseta(),electron->pT());
	return (teff);      

}



double MuonEff(HEP::Particle* muon) {
        static HEP::BinnedFn2D<double> _muEff2d({{0,0.1,2.5,2.7,DBL_MAX}}, //< |eta|
                                                {{0,5.0,10.0,DBL_MAX}}, //< pT
                                                {{0.0,0.9,0.9, // |eta| < 0.1
                                                  0.0,0.98,0.99, // |eta| < 2.5
                                                  0.0,0.0,0.99, // |eta| < 2.7
                                                  0.0,0.0,0.0 // large eta
                                                }});


          double teff=_muEff2d.get_at(muon->abseta(),muon->pT());
          return (teff);  
        }




HEP::P4 smearMET(const HEP::P4 &pTmiss, std::mt19937 &engine) {

  // uncertainty on MET is ~ 20 GeV on perp and , assume normal distribution
  double MET=pTmiss.pT();
  std::normal_distribution<> d(MET,20.0);

  double smear_factor=d(engine);
  if(smear_factor <=0.0) return pTmiss;

  return HEP::P4::mkEtaPhiMPt(pTmiss.eta(),pTmiss.phi(),pTmiss.m(), smear_factor);

}


*/
