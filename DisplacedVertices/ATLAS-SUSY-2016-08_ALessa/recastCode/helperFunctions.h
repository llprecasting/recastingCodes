//Some helper functions for computing relevant quantities
#include "dataReader.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


using namespace Pythia8;


class DisplacedVertex: public Particle {
  public:

    DisplacedVertex(){}
    DisplacedVertex(const Particle& pt) : Particle(pt) {}
    virtual ~DisplacedVertex() {}

    vector<Particle*> decayProducts;
    double DVeff;
	int nTracks() const {return decayProducts.size();}
	double mDV() const {
        double px = 0.0, py = 0.0, pz = 0.0, e = 0.0;
		for (int i=0; i < decayProducts.size(); ++ i){
			px += decayProducts[i]->px();
            py += decayProducts[i]->py();
            pz += decayProducts[i]->pz();
            // e += decayProducts[i]->e();
            //Use pion mass assumption:
            e += sqrt(decayProducts[i]->pAbs2() + pow(0.139,2));
		}
		Vec4 pCharged(px,py,pz,e);
		return pCharged.mCalc();
    }

 };

double d0Calc(Vec4 xv, Vec4 p)
//Compute the transverse impact parameter for a particle with momentum p created at vertex xv
//assuming no magnetic field
{
    double calpha = p.px()/p.pT();
    double salpha = p.py()/p.pT();

    double a0 = xv.px()*salpha - xv.py()*calpha;

    return fabs(a0);
}

double getEvEff(double MET, double Rmax)
//Compute the event efficiency using information
//from ATLAS efficiency maps
{

    vector< vector<double> > effData;
    string table;
    if (Rmax < 1150.){
        effData = effEv_table22; //Fig.20a
    }
    else if (Rmax < 3870.){
        effData = effEv_table23; //Fig.20b
    }
    else{effData = effEv_table24;} //Fig.20c

    double METbinLow,METbinHigh,eff;
    for (int i = 0; i < effData.size(); ++i){
        METbinLow = effData[i][1];
        METbinHigh = effData[i][2];
        eff = effData[i][3];
        if (METbinLow <= MET && MET < METbinHigh){return eff;}
    }

    return 0.;
}

double getDVEff(double mDV, int nTracks, double Rmax)
//Compute the DV reconstruction efficiency using information
//from ATLAS efficiency maps
{

    vector< vector<double> > effData;
    string table;
    if ( Rmax < 4.){return 0.;}
    else if (Rmax < 22.){
        effData = effDV_table25; //Fig.21a
    }
    else if (Rmax < 25.){
        effData = effDV_table26; //Fig.21b
    }
    else if (Rmax < 29.){
        effData = effDV_table27; //Fig.21c
    }
    else if (Rmax < 38.){
        effData = effDV_table28; //Fig.21d
    }
    else if (Rmax < 46.){
        effData = effDV_table29; //Fig.21e
    }
    else if (Rmax < 73.){
        effData = effDV_table30; //Fig.21f
    }
    else if (Rmax < 84.){
        effData = effDV_table31; //Fig.22a
    }
    else if (Rmax < 111.){
        effData = effDV_table32; //Fig.22b
    }
    else if (Rmax < 120.){
        effData = effDV_table33; //Fig.22c
    }
    else if (Rmax < 145.){
        effData = effDV_table34; //Fig.22d
    }
    else if (Rmax < 180.){
        effData = effDV_table35; //Fig.22e
    }
    else if (Rmax < 300.){
        effData = effDV_table36; //Fig.22f
    }
    else{return 0.;}


    double mDVbinLow,mDVbinHigh,nTracksbinLow,nTracksbinHigh,eff;
    for (int i = 0; i < effData.size(); ++i){
        mDVbinLow = effData[i][1];
        mDVbinHigh = effData[i][2];
        nTracksbinLow = int(effData[i][4]);
        nTracksbinHigh = int(effData[i][5]);
        eff = effData[i][6];
        if (mDVbinLow <= mDV && mDV < mDVbinHigh &&
                nTracksbinLow <= nTracks && nTracks < nTracksbinHigh){
            return eff;}
    }

    return 0.;
}

std::vector<Particle*> getDaughters(Particle* &mom, Event &event)
//Returns a vector with all the (stable) daughters of a given particle
{

    //Collect daughter indices:
    std::vector<int> daughtersAll = mom->daughterListRecursive();
    std::vector<int> dIndices;
    for (int j = 0; j < daughtersAll.size(); ++j){
        if (!event[daughtersAll[j]].isFinal()) continue;
        if (std::find(dIndices.begin(),dIndices.end(),daughtersAll[j]) != dIndices.end()) continue;
        dIndices.push_back(daughtersAll[j]);
    }
    //Get list of daughters:
    std::vector<Particle*> daughters;
    for (int j = 0; j < dIndices.size(); ++j){
        daughters.push_back(&event[dIndices[j]]);
    }

    return daughters;
}

std::vector<Particle*> getRhadrons(Event &event)
//Returns a vector with the Rhadrons in the event
{

    std::vector<Particle*> Rhadrons;
    Rhadrons.resize(0);
    // Loop over final R-hadrons in the event: kinematic distribution
    for (int i = 0; i < event.size(); ++i) {
      //R-hadrons
      if (abs(event[i].status()) == 104) {
        Rhadrons.push_back(&event[i]);
      }
    }

    return Rhadrons;
}


std::vector<DisplacedVertex> getDVs(Event &event,double minTrackPT, double minTrackD0){
//Returns a vector with displaced vertices satisfying minimum selection criterium

    std::vector<DisplacedVertex> DVs;  //Store DVs that passes the cuts
    //Get Rhadrons from event:
    std::vector<Particle*> Rhadrons = getRhadrons(event);
    int nRhadrons = Rhadrons.size();
    for (int i = 0; i < nRhadrons; ++i) {
        Vec4 RhadronVertex = Rhadrons[i]->vDec();
        //Get daughters from R-hadron:
        std::vector<Particle*> daughters = getDaughters(Rhadrons[i],event);
        //Define DV candidate object to store candidate tracks as well:
        DisplacedVertex DV(*Rhadrons[i]);
        //Get selected decay products:
        Vec4 pCTotal;
        Vec4 pC;
        int nGoodTracks = 0;
        double d0,z0,DVeff,Rdecay;
        for (int j=0; j < daughters.size(); ++j){
            if(daughters[j]->isNeutral()) continue; //Only consider charged particles
            pC = daughters[j]->p();
            if (pC.pT() < minTrackPT) continue; //p_T > 1 GeV for tracks
            d0 = d0Calc(RhadronVertex,pC);
            if (fabs(d0) < minTrackD0) continue;  //|d_0| > 2mm for tracks
            DV.decayProducts.push_back(daughters[j]);
        }
        DVs.push_back(DV);
    } //End of loop over R-hadrons
    return DVs;
}


Vec4 getMissingMomentum(Event &event)
//Returns the missing 4-vector for the event (sums up neutralinos and neutrinos)
{
    Vec4 missingETvec;
    for (int i = 0; i < event.size(); ++i) {
      // Final state only
      if (!event[i].isFinal()) continue; //Ignore intermediate states
      if (event[i].isVisible()) continue; //Ignore intermediate states
      if (!event[i].isNeutral()) continue; //Ignore charged particles
      if (event[i].isHadron()) continue;  //Ignore usual hadrons
      if (event[i].colType() != 0) continue; //Ignore charged particles
      if (event[i].status() == 104) continue; //Ignore stable R-Hadrons
      if (event[i].idAbs() == 22) continue; //Ignore photons
      // Missing momentum
      missingETvec += event[i].p();
    }
    return missingETvec;
}




bool applyJetCuts(Event &event, fastjet::JetDefinition jetDef,
        double pTjet, double maxJetChargedPT,
        double minJetPt1, double minJetPt2,    double minPVdistance)
//Apply basic selection cuts based on jet requirements
{
    //Get jets:
    vector <fastjet::PseudoJet> inclusiveJets, jets;
    std::vector <fastjet::PseudoJet> fjInputs; //particles for applying jet clustering
    // Loop over particles in the event: kinematic distribution
    for (int i = 0; i < event.size(); ++i) {
        // Require visible particles inside detector.
        if (!event[i].isFinal()){continue;}
        fastjet::PseudoJet particleTemp(event[i].px(),event[i].py(),event[i].pz(),event[i].e());
        particleTemp.set_user_index(i);
        fjInputs.push_back(particleTemp);
    }
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(pTjet);
    jets    = sorted_by_pt(inclusiveJets);

    int nGoodJets = 0;
    bool passJetCuts = false;
    double pTCharged;
    Particle jetParticle;
    for(int i = 0; i < jets.size(); ++i){
        pTCharged = 0.;
        vector<fastjet::PseudoJet> constituents = jets[i].constituents();
        for (int j = 0; j < constituents.size(); j++) {
            jetParticle = event[constituents[j].user_index()];
            if (jetParticle.isNeutral()){continue;}
            if (jetParticle.vProd().pT() > minPVdistance){continue;}
            pTCharged += jetParticle.pT();
        }
        if (pTCharged > maxJetChargedPT) {continue;}
        if (jets[i].pt() > minJetPt1){passJetCuts = true; break;}
        if (jets[i].pt() > minJetPt2){++nGoodJets;}

    }
    if (nGoodJets >= 2){passJetCuts = true;}
    if (!passJetCuts){return false;}

    return true;

}
