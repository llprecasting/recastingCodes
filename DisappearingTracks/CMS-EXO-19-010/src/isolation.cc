#include "include/isolation.h"

double d0Calc(const HEP::P4 &xv, const HEP::P4 &p)
//Compute the transverse impact parameter for a particle with momentum p created at vertex xv
//assuming no magnetic field
// Look at x-y plane and find clsost distance: this is |xv| sin varphi, where the angle is between the vertex
// and the momentum
// this is just given by xv cross phat, where each is represented as a 2d vector
{
    double calpha = p.px()/p.pT();
    double salpha = p.py()/p.pT();

    double a0 = xv.px()*salpha - xv.py()*calpha;

    return fabs(a0);
}

double dzCalc(const HEP::P4 &xv, const HEP::P4 &p)
//Compute the longitudinal impact parameter for a particle with momentum p created at vertex xv
//assuming no magnetic field
// dz is the z value at the point of closest d0.
// v_min = xv - lambda \hat{p}
// lambda = |xv| cos varphi = xv dot \hat{p} (in 2d)
{
    //double calpha = p.px()/p.pT();
    //double salpha = p.py()/p.pT();

    //double a0 = xv.px()*salpha - xv.py()*calpha;
    double lambda= fabs((p.px()*xv.px()+p.py()*xv.py()))/p.pT();
    double zmin = xv.pz()-lambda*p.pz()/p.pT();
    return fabs(zmin);
}

double dz_sintheta(const HEP::P4 &xv, const HEP::P4 &p)
//Compute the longitudinal impact parameter for a particle with momentum p created at vertex xv
{
    //double calpha = p.px()/p.pT();
    //double salpha = p.py()/p.pT();

    //double a0 = xv.px()*salpha - xv.py()*calpha;
    if(p.p() == 0.0) return 0.0;
    double dz = dzCalc(xv,p);
    
    return fabs(dz*p.pT()/p.p());
    //return fabs(xv.p());
 /*
 if(xv.p() == 0.0) return 0.0;
    double dz = dzCalc(xv,p);

    return fabs(dz*xv.pT()/xv.p());
    */
   
}



//template<typename T1, typename T2> std::vector<T1*> overlapRemoval

//isolation 

double sumpTisolation(const HEP::P4 &p4, const HEP::Event* evt, const double &DR)
{
    // compute pT around a cone of size DR around p4, for charged tracks -> charged hadrons, electrons, muons
    
    
    double sumpT=0.0;
    
    for(auto chad : evt->get_charged_hadrons())
    {
      //if(chad->does_decay()) continue; // this is now fixed upstream
        if(p4.deltaR_eta(chad->mom()) > DR) continue;
        double tPt=chad->pT();
        if(tPt > 0.1)
        {
            sumpT+=tPt;
        }
    }

    for(auto hscp : evt->HSCPs())
    {
        if(p4.deltaR_eta(hscp->mom()) > DR) continue;
        double tPt=hscp->pT();
        if(tPt > 0.1)
        {
            sumpT+=tPt;
        }
    }


    for(auto elec : evt->electrons())
    {
        if(p4.deltaR_eta(elec->mom()) > DR) continue;
        
        double tPt=elec->mom().pT();
        if(tPt > 0.1)
        {
            sumpT+=tPt;
        }
    }

    for(auto muon : evt->muons())
    {
        if(p4.deltaR_eta(muon->mom()) > DR) continue;

        double tPt=muon->mom().pT();
        if(tPt > 0.1)
        {
            sumpT+=tPt;
        }
    }
    // Generally we should subtracting itself, 
    // to be general (e.g. if we want to do this for photon isolation) we could skip that step but for ease here I retain it

    sumpT-=p4.pT();
    return sumpT;


}


double sumETIsolation(const HEP::P4 &p4, const HEP::Event* evt, const double &DR)
{
    // compute ET around a cone of size DR around p4, for everything -> charged hadrons, neutral hadrons, electrons, muons, photons, subtracting itself
    double sumET =0.0;
    
    for(auto chad : evt->get_charged_hadrons())
    { 
      //if(chad->does_decay()) continue;
        if(p4.deltaR_eta(chad->mom()) > DR) continue;

        sumET+=HEP::ET(*chad);
    }

    for(auto neut : evt->get_neutral_hadrons())
    {
        if(p4.deltaR_eta(*neut) > DR) continue;

        sumET+=HEP::ET(*neut);
    }

    for(auto hscp : evt->HSCPs())
    {
        if(p4.deltaR_eta(hscp->mom()) > DR) continue;

        sumET+=HEP::ET(hscp->mom());

    }

    for(auto elec : evt->electrons())
    {
        if(p4.deltaR_eta(elec->mom()) > DR) continue;

        sumET+=HEP::ET(elec->mom());
    }

    for(auto muon : evt->muons())
    {
        if(p4.deltaR_eta(muon->mom()) > DR) continue;

        sumET+=HEP::ET(muon->mom());
    }

     for(auto photon : evt->photons())
    {
        if(p4.deltaR_eta(photon->mom()) > DR) continue;

        sumET+=HEP::ET(photon->mom());
    }

    sumET-=HEP::ET(p4);

    return sumET;


};




double getcos(const double eta)
{
  if(eta ==0.) return 0.0;
  double myexp=exp(-2.0*eta);
  return (1.0-myexp)/(1.0+myexp);

}


double coneAreaApprox(const HEP::P4 &p4, const double &DR)
{
  // Gives an approximation for the area of a cone within fixed DR
  double eta0=p4.abseta();
  double rectarea=2.0* DR*(2.0/(1.0 + exp(2.0*(eta0 - DR))) - 2.0/(1.0 + exp(2.0*(eta0 + DR))));
  double res = M_PI/4.0*rectarea;

  return res;
}



double sumCaloE(const HEP::P4 &p4, const HEP::Event* evt, const double &DR)
{
    // compute All calorimeter energy around a cone of size DR around p4, for everything -> charged hadrons, neutral hadrons, electrons, muons, photons, subtracting itself
    double sumET =0.0;
    double calE;
    double threshold=0.5;
    for(auto chad : evt->get_charged_hadrons())
    {
      //if(chad->does_decay()) continue;
        if(p4.deltaR_eta(chad->mom()) > DR) continue;

        // check actually registers in HCAL
        if(chad->does_decay())
        {
          //abs(ptc.vProd().pz()) > 5500. || ptc.vProd().pT() > 2950.
          //if((chad->decay_vertex().pT() < 2950.)  && (abs(chad->decay_vertex().pz()) < 5500.)) continue;         
          if((chad->decay_vertex().pT() < 1100.)  && (abs(chad->decay_vertex().pz()) < 2700.)) continue;        
        }
        else
        {
          if((chad->prod_vertex().pT() > 2950.)  || (abs(chad->prod_vertex().pz()) > 5500.)) continue; 
        }

        calE=chad->E();
       

        if(calE > threshold)
        {
          //std::cout << "Charged hadron nearby!" << calE << " " << chad->pid() << std::endl;
          sumET+=calE;

           if(calE > 10)
         {
            std::cout << "Energetic charged hadron " << chad->pid() << ", " << calE << std::endl;
          }
        }
    }

    for(auto neut : evt->get_neutral_hadrons())
    {
        if(p4.deltaR_eta(*neut) > DR) continue;

        calE=neut->E();


        if(calE > threshold)
        {
          sumET+=calE;
           if(calE > 10)
        {
          std::cout << "Energetic neutral hadron " << calE << std::endl;
        } 

        }
    }
    
    /*
    for(auto hscp : evt->HSCPs())
    {
        if(p4.deltaR_eta(hscp->mom()) > DR) continue;
         if(hscp->does_decay())
        {
          //abs(ptc.vProd().pz()) > 5500. || ptc.vProd().pT() > 2950.
          //if((chad->decay_vertex().pT() < 2950.)  && (abs(chad->decay_vertex().pz()) < 5500.)) continue;         
          if((hscp->decay_vertex().pT() < 1100.)  && (abs(hscp->decay_vertex().pz()) < 2700.)) continue;        
        }
        else
        {
          if((hscp->prod_vertex().pT() > 2950.)  || (abs(hscp->prod_vertex().pz()) > 5500.)) continue; 
        }


        calE=hscp->E();
        if(calE > threshold)
        {
          sumET+=calE;
        }

    }
    */
    for(auto elec : evt->electrons())
    {
        if(p4.deltaR_eta(elec->mom()) > DR) continue;

        calE=elec->E();
        if(calE > threshold)
        {
          sumET+=calE;
        }
    }

    /*  Muons deposit very little energy
    for(auto muon : evt->muons())
    {
        if(p4.deltaR_eta(muon->mom()) > DR) continue;

        calE=muon->E();
        if(calE > threshold)
        {
          sumET+=calE;
           if(calE > 30)
        {
          std::cout << "Energetic muon " << calE << std::endl;
        } 
        }
    }

    */
     for(auto photon : evt->photons())
    {
        if(p4.deltaR_eta(photon->mom()) > DR) continue;

        calE=photon->E();
        if(calE > threshold)
        {
          sumET+=calE;
        }
    }
      
     //sumET-=p4.E();
     return sumET;


};




bool IsLooseTightLepton(const  HEP::Particle* Lep, const HEP::Event* event,
		   const double DeltaRp, const double DeltaRE, const double pratio, const double Eratio)
{
  return true;
  // varconeXX for momentum, coneXX for transverse energy
  double pt = Lep->pT();
  if(pt < 1.0) // just to be sure ... cut is much stricter than this on pt anyway
    {
      return false;
    };
  double DRp = std::min(10.0/pt,DeltaRp);
  //double imini = (PHYSICS->Isol->eflow->sumIsolation(Lep,event.rec(),DRp,1.0,IsolationEFlow::ALL_COMPONENTS))/pt;
  
    double imini = sumpTisolation(Lep->mom(), event, DRp)/pt;

  if (imini > pratio)
    {
      return false;

    };
  
  double jmini = sumETIsolation(Lep->mom(),event,DeltaRE)/pt;
  
  if (jmini > Eratio)
    {
      return false;

    };


  return true;
}


bool IsHighPTCaloOnly(const HEP::Particle* Lep, const HEP::Event* event,
		    const double DeltaRE)
{
  return true;
  // double eta = std::fabs(Lep->eta());
  // double pt = Lep->pt();

  // double DR = 10./std::min(std::max(pt,50.),200.);
  // double imini = (PHYSICS->Isol->eflow->sumIsolation(Lep,
  //   event.rec(),DR,0.,IsolationEFlow::ALL_COMPONENTS))/pt;

  // if( eta<etaTH && pt>ptTH && imini<0.4) return true;

  // varconeXX for momentum, coneXX for transverse energy
  double pt = Lep->pT();
  if(pt < 1.0) // We only apply this for pt> 200, but here we're just being safe
    {
      return false;
    };
  //double DRp = std::min(10.0/pt,DeltaRp);
  double ETratio=std::max(0.015*pt,3.5);
  
  double jmini = sumETIsolation(Lep->mom(),event,DeltaRE);
 
  if (jmini > ETratio)
    {
      return false;

    };


  return true;
}








