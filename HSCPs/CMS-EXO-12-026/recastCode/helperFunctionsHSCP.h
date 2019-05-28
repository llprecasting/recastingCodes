//Some helper functions for computing relevant quantities
#include <string>
#include <vector>
#include <vector>
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TH3.h"


using namespace Pythia8;


std::vector<Particle*> getHSCPs(Event &event)
//Returns a vector with the HSCPs in the event
{

    std::vector<Particle*> HSCPs;
    HSCPs.resize(0);
    // Loop over final particles in the event
    for (int i = 0; i < event.size(); ++i) {
      //HSCPs
      if (!event[i].isFinal()) continue;
      if (abs(event[i].m()) < 10.) continue;
      if (abs(event[i].status()) == 104 || event[i].isCharged()) {
        HSCPs.push_back(&event[i]);
      }
    }
    return HSCPs;
}



bool isIsolated(Event &event, Particle* HSCP)
{

    double caloIso=0.0;
    double tkIso=0.0;
    double deltaR;
    for(int i=0; i < event.size() ; i++){
        if (!event[i].isFinal()) continue;
        int absPdg = abs(event[i].id());
        if (i == HSCP->index()) continue; //skip the candidate
        if(absPdg==12 || absPdg==14 || absPdg==16) continue; //skip neutrinos
        deltaR = sqrt(pow(event[i].eta()-HSCP->eta(),2) + pow(event[i].phi()-HSCP->phi(),2));
        if(deltaR > 0.3) continue;
        if(event[i].isCharged()) tkIso += event[i].pT();
        caloIso += event[i].e();
    }
    if (caloIso/HSCP->pAbs() > 0.3)return false;
    if (tkIso > 50.) return false;
    return true;
}

std::vector<Particle*> getIsolatedHSCPs(Event &event){
//Returns a vector with displaced vertices satisfying minimum selection criterium

    std::vector<Particle*> HSCPs = getHSCPs(event);
    std::vector<Particle*> IsoHSCPs;
    int nHSCPs = HSCPs.size();
	for (int i = 0; i < nHSCPs; ++i) {
        if (isIsolated(event,HSCPs[i])) IsoHSCPs.push_back(HSCPs[i]);
    }
    return IsoHSCPs;
}


double computeFlong(Particle* HSCP, double width = 0.){
// Compute the fraction of long-lived particles for a given width

    double ctau = -1.0;
    if (width > 0.) ctau = 1.975e-16/width;  //c*tau in meters
    double x = 11.0; // CMS detector size
    if(fabs(HSCP->eta())<0.8) x = 9.0;
    else if(fabs(HSCP->eta())<1.1) x = 10.0;
    double Flong = ctau>0?exp(-HSCP->m0()*x*(1.0/ctau)/HSCP->pAbs()):1.0; //Compute fraction of long-lived
    
    return Flong;
}

std::vector< std::pair<double,double> > computeHSCPeff(Particle* HSCP,TH3F* histoTrigger,
											TH3F** histoOffline){
   //Compute the HSCP efficiency for each mass cut

    int nM = 4;

    //Total efficiencies
    std::vector<double> totalEff(nM,0.);
    std::vector<double> totalEffErr(nM,0.);
    std::vector< std::pair<double,double> > result(nM+1,std::make_pair(0.,0.));

    //Trigger probabilities and error:
    double ProbTrigger = 0., ProbTriggerErr = 0.;
    //Online probabilities and error for each signal region/mass cut:
    std::vector<double> ProbOnline(nM,0.),ProbOnlineErr(nM,0.);

    if(HSCP->pT() < 40. || fabs(HSCP->eta()) > 2.1) return result;//don't waste time on stuff that will never make it
    int AbsPdg = abs(HSCP->id());
    if((AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 
        || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 
        || AbsPdg==1093214 || AbsPdg==1093324)) return result; //Skip neutral gluino RHadrons
    if((AbsPdg==1000622 || AbsPdg==1000642 
        || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 
        || AbsPdg==1006333)) return result;  //skip neutral stop RHadrons

    int BinX = histoTrigger->GetXaxis()->FindBin(HSCP->pT());
    int BinY = histoTrigger->GetYaxis()->FindBin(HSCP->pAbs()/HSCP->e()); 
    int BinZ = histoTrigger->GetZaxis()->FindBin(fabs(HSCP->eta()));
    ProbTrigger = histoTrigger->GetBinContent(BinX, BinY, BinZ);
    ProbTriggerErr = histoTrigger->GetBinError(BinX, BinY, BinZ);
    for(unsigned int Mi=0;Mi<nM;Mi++){
        //Only use efficiencies for mHSCP*0.6 > mCUT (mCUT = 100*Mi)
        if ((HSCP->m0())*0.6 > 100.*Mi){
        	ProbOnline[Mi] = histoOffline[Mi]->GetBinContent(BinX, BinY, BinZ);
        	ProbOnlineErr[Mi] = histoOffline[Mi]->GetBinError(BinX, BinY, BinZ);
        }
    }

    result[0] = std::make_pair(ProbTrigger,ProbTriggerErr);
    for(unsigned int Mi=0;Mi<nM;Mi++){
        result[Mi+1] = std::make_pair(ProbOnline[Mi],ProbOnlineErr[Mi]);
    }

    return result;
}

std::vector< std::pair<double,double> > computeTotalEff(std::vector< std::vector< std::pair<double,double> > > effs,
															std::vector<double> Flong){
   //Given the efficiencies for the HSCPs and the fraction of long-lived particles
//for each HSCP compute the total event efficiency


    int nM = 4;
    //Total efficiencies
    std::vector<double> totalEff(nM,0.);
    std::vector<double> totalEffErr(nM,0.);
    std::vector< std::pair<double,double> > result(nM,std::make_pair(0.,0.));

    //Number of isolated HSCPs:
    int nhscp = effs.size();
    if (nhscp < 1) return result;
    if (nhscp > 2){
        cout << "Can not handle more than 2 HSCPs per event yet" << endl;
        return result;
    }


    //Trigger probabilities, error and fraction of long-lived particles:
    std::vector<double> ProbTrigger(2,0.),ProbTriggerErr(2,0.),FlongV(2,0.);
    //Online probabilities and error for each signal region/mass cut:
    std::vector< std::vector<double> > ProbOnline,ProbOnlineErr;
    ProbOnline.resize(2, std::vector<double>(nM, 0.));
    ProbOnlineErr.resize(2, std::vector<double>(nM, 0.));

    for (int ihscp=0; ihscp < effs.size(); ++ihscp){
    	FlongV[ihscp] = Flong[ihscp];
        ProbTrigger[ihscp] = effs[ihscp][0].first;
        ProbTriggerErr[ihscp] = effs[ihscp][0].second;
        for(unsigned int Mi=0;Mi<nM;Mi++){            
            ProbOnline[ihscp][Mi] = effs[ihscp][Mi+1].first;
            ProbOnlineErr[ihscp][Mi] = effs[ihscp][Mi+1].second;
        }
    }

    double totalTrigger = 0.;
    double totalTriggerErr = 0.;
    //Compute total probabilities for both HSCPs (if there is only one it is equivalent to using its prob.)
    totalTrigger = ProbTrigger[0] + ProbTrigger[1] - ProbTrigger[0]*ProbTrigger[1];
    totalTriggerErr = pow(ProbTriggerErr[0]*(1.-ProbTrigger[1]),2) + pow(ProbTriggerErr[1]*(1.-ProbTrigger[0]),2);
    for(unsigned int Mi=0;Mi<nM;Mi++){
        double total = totalTrigger;
        double totalErr = totalTriggerErr;
        double FlongEff = 0.;
        total *= (ProbOnline[0][Mi] + ProbOnline[1][Mi] - ProbOnline[0][Mi]*ProbOnline[1][Mi]);
        if (total != 0.){
			FlongEff = FlongV[0]*FlongV[1];
			FlongEff += FlongV[0]*(1.-FlongV[1])*ProbTrigger[0]*ProbOnline[0][Mi]/(total);
			FlongEff += FlongV[1]*(1.-FlongV[0])*ProbTrigger[1]*ProbOnline[1][Mi]/(total);
        }
        totalErr += pow(ProbOnlineErr[0][Mi]*(1-ProbOnline[1][Mi]),2) + pow(ProbOnlineErr[1][Mi]*(1-ProbOnline[0][Mi]),2);
        totalEff[Mi] = total*FlongEff;
        totalEffErr[Mi] = totalErr;        
        result[Mi] = std::make_pair(totalEff[Mi],totalEffErr[Mi]);
     }   
    
    return result;
}



