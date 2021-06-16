

#include "HSCP_ATLAS.h"
#include "include/isolation.h"

const std::vector<std::string> allSRs={"CAND1_0","CAND1_1","CAND1_2","CAND1_3","CAND2_0","CAND2_1","CAND2_2","CAND2_3"};
const std::vector<std::string> CAND1SRs={"CAND1_0","CAND1_1","CAND1_2","CAND1_3"};
const std::vector<std::string> CAND2SRs={"CAND2_0","CAND2_1","CAND2_2","CAND2_3"};

const std::vector<std::string> CAND1mTOFS={"mTOF>175","mTOF>375","mTOF>600","mTOF>825"};
const std::vector<std::string> CAND2mTOFS={"mTOF>150","mTOF>350","mTOF>575","mTOF>800"};

std::vector<double> massToF2Cand_min  = {150.,350.,575.,800.};
std::vector<double> massToF1Cand_min  = {175.,375.,600.,825.};

void HSCP_ATLAS::init() {
   // cout << "initialising analysis" << endl;


  cout << "---------------------------------------------------" <<std::endl;
  cout << "-- ATLAS HSCP search 36.1 fb^-1                  --" <<std::endl;
  cout << "-- arXiv:1902.01636                              --" << std::endl;
  cout << "-- By Mark Goodsell (goodsell@lpthe.jussieu.fr)  --" << std::endl;
  cout << "-- Based on HEPData material                     --" << std::endl;
  cout << "-- and LLPRecastingCodes github code by A. Lessa --" << std::endl;
  cout << "---------------------------------------------------" <<std::endl;



for(std::string region : allSRs)
      {
	AddRegionSelection(region);
      }

AddCut("Trigger",allSRs);
//AddCut("NTracks > 1",allSRs);
AddCut("Preselection",allSRs);
//AddCut("|eta| < 2",allSRs);

AddCut("2Loose",CAND2SRs);
AddCut("1Tight",CAND1SRs);

// std::vector<float> massToF2Cand_min  = {150.,350.,575.,800.};
  
AddCut("mTOF>150","CAND2_0");
AddCut("mTOF>350","CAND2_1");
AddCut("mTOF>575","CAND2_2");
AddCut("mTOF>800","CAND2_3");

// std::vector<float> massToF1Cand_min  = {175.,375.,600.,825.};
AddCut("mTOF>175","CAND1_0");
AddCut("mTOF>375","CAND1_1");
AddCut("mTOF>600","CAND1_2");
AddCut("mTOF>825","CAND1_3");

    


/*
      //Load the ATLAS histograms:
  TFile* InputFile = new TFile("/home/mark/Data/SModelS/EWStudy/HSCP/200331/MyHSCP/recastCode/ATLAS_data/HEPData-ins1718558-v2-root.root");
  TH1F* EtmissTurnOn = (TH1F*) InputFile->GetDirectory("Table 22")->Get("Hist1D_y1");
  TH2F* SingleMuTurnOn = (TH2F*) InputFile->GetDirectory("Table 23")->Get("Hist2D_y1");
  TH2F* LooseEff = (TH2F*) InputFile->GetDirectory("Table 25")->Get("Hist2D_y1");
  TH2F* TightPromotionEff = (TH2F*) InputFile->GetDirectory("Table 26")->Get("Hist2D_y1");
  TH1F* MToFFullDet = (TH1F*) InputFile->GetDirectory("Table 29")->Get("Hist1D_y1");
  TH1F* MToFFullDetErr = (TH1F*) InputFile->GetDirectory("Table 29")->Get("Hist1D_y1_e1");
*/
/*
   Efficiency1D EtmissTurnOn("HSCP-ATLAS/Table22.csv");   
   Efficiency2D SingleMuTurnOn("HSCP-ATLAS/Table23.csv");
   Efficiency2D LooseEff("HSCP-ATLAS/Table25.csv");   
   Efficiency2D TightPromotionEff("HSCP-ATLAS/Table26.csv");  

   Efficiency1D MToFFullDet("HSCP-ATLAS/Table29.csv");
   Efficiency1D MToFFullDetErr("HSCP-ATLAS/Table29.csv",true); // second argument means return the uncertainty
      */

   EtmissTurnOn.readcsv("HSCP-ATLAS/Table22.csv");   
   SingleMuTurnOn.readcsv("HSCP-ATLAS/Table23.csv");
   LooseEff.readcsv("HSCP-ATLAS/Table25.csv");   
   TightPromotionEff.readcsv("HSCP-ATLAS/Table26.csv");  

   MToFFullDet.readcsv("HSCP-ATLAS/Table29.csv");
   MToFFullDetErr.readcsv("HSCP-ATLAS/Table29.csv",true); // second argument means return the uncertainty
  

}


HSCP_ATLAS::~HSCP_ATLAS() { 
  

};


bool decayInsideAtlas(HEP::Particle * p)
{

      if(!(p->does_decay())) return false;
      //cout << "Doens't decay!" << endl;
      //cout << p->decay_vertex().pT()  << ", " << fabs(p->decay_vertex().pz()) << endl;
      if((p->decay_vertex().pT() < 12e3) && (fabs(p->decay_vertex().pz()) < 23e3)) return true;

      //cout << "is fine " << endl;
      return false;

}


void HSCP_ATLAS::Execute(std::mt19937 &engine) {

      
      static std::uniform_real_distribution<double> rd(0.0,1.0);
      
      double Etmiss=Event->met();
      bool TriggerAccept = false;
      //cout << "MET: " << Etmiss << "  " ;
      if (Etmiss > 300.) 
      {
            TriggerAccept = true;
      }
      else
      {

      double eff_Met = EtmissTurnOn.get_at(Etmiss);
	double metRandom = rd(engine);
	if (metRandom < eff_Met) TriggerAccept = true;
      }

      std::vector<HEP::Particle*> initial_hscps;

     // cout << Event->HSCPs().size() << " initial HSCPs" << endl;
      
            for(auto particle : Event->HSCPs())
            {
                  
                  if (decayInsideAtlas(particle)) continue;
                  //Only single charged particles are considered
                  if (abs(particle->get_3Q()) != 3) continue;
                  
                  initial_hscps.push_back(particle);
                  //If the MET trigger was not activated, check for the muon trigger:
                  if (!TriggerAccept)
                         {
                        double abseta  = particle->abseta();
			      double beta  = particle->mom().p()/particle->E();
			
			      double effTrig   = SingleMuTurnOn.get_at(abseta, beta);
			      double triggerRandom = rd(engine);
			      if (triggerRandom < effTrig) {
				TriggerAccept = true;
				
			      }
                  }    
            }
      
      //cout << initial_hscps.size() << " passed HSCPs" << endl;
      ApplyCut(TriggerAccept,"Trigger");
     
      if(!TriggerAccept) return;
      // bool ntracks = initial_hscps.size();
       //ApplyCut(ntracks>)

      std::vector<HEP::Particle*> looseCandidates;
	std::vector<HEP::Particle*> tightCandidates;

      int pass_preselection=0;
      for(auto particle : initial_hscps)
      {
            double abseta  = particle->abseta();
            double p = particle->mom().p();
            double pT = particle->pT();
		double beta  = p/particle->E();
            if (abseta > 2.0) continue;
		if (pT < 70.0) continue;
		if (p < 100.0) continue;

            pass_preselection++;
            if(beta > 0.95) continue;
            // Estimate efficiencies
		double effLoose = LooseEff.get_at(abseta, beta);
		double effPromotion = TightPromotionEff.get_at(abseta, beta);
		double looseRandom = rd(engine);
		double promoteRandom = rd(engine);

		//avgCandEff += effLoose;

		if (looseRandom < effLoose) {
		    looseCandidates.push_back(particle);

                // Tight selection
                if(abseta > 1.65) continue;
                if(beta > 0.8) continue;
		    if (particle->mom().p() < 200.) continue;
		    // Sample tight promotion of candidate
		    if (promoteRandom < effPromotion) {
		    	tightCandidates.push_back(particle);
		    }
		  }
      }

      ApplyCut(pass_preselection>0,"Preselection");
      //cout << pass_preselection << " passed p, then " << looseCandidates.size() << " and  " << tightCandidates.size() << endl;
      ApplyCut(looseCandidates.size() == 2,"2Loose");
      ApplyCut(tightCandidates.size() == 1,"1Tight");
      //std::vector<int> nEvts_SR1Cand = {0,0,0,0};
      //std::vector<int> nEvts_SR2Cand = {0,0,0,0};

      
      if(looseCandidates.size() == 2)
      {
            double mass1 = fabs(looseCandidates[0]->mass());
		double mass2 = fabs(looseCandidates[1]->mass());
	    // Events with exactly two loose candidates
	    // Sample the ToF mass for two full-detector candidates
	    
          
	      double massToF_mean1       = MToFFullDet.get_at(mass1);
            double massToF_resolution1 = MToFFullDetErr.get_at(mass1);

            double massToF_mean2       = MToFFullDet.get_at(mass2);
            double massToF_resolution2 = MToFFullDetErr.get_at(mass2);

	    std::normal_distribution<double> gaussToF1(massToF_mean1, massToF_resolution1);
	    std::normal_distribution<double> gaussToF2(massToF_mean2, massToF_resolution2);
	    double massToF1   = gaussToF1(engine);
	    double massToF2   = gaussToF2(engine);

		for (int j = 0; j < massToF2Cand_min.size(); ++j){
			// Apply final mass requirements for each SR
                  bool pass_cut=true;
			if (std::min(massToF1, massToF2)  < massToF2Cand_min[j]) pass_cut=false;

                  ApplyCut(pass_cut,CAND2mTOFS[j]);
			
		}
      
            
      }
      else  if(tightCandidates.size() == 1)
      {
            double mass1 = fabs(tightCandidates[0]->mass());
		
	      double massToF_mean1       = MToFFullDet.get_at(mass1);
            double massToF_resolution1 = MToFFullDetErr.get_at(mass1);

	    std::normal_distribution<double> gaussToF1(massToF_mean1, massToF_resolution1);

	    double massToF1   = gaussToF1(engine);
           // cout << "massToF1: " << massToF1 << ": ";

		for (int j = 0; j < massToF1Cand_min.size(); ++j){
			// Apply final mass requirements for each SR
                  //cout << massToF1Cand_min[j] << ", ";
                  bool pass_cut=true;
			if (massToF1  < massToF1Cand_min[j]) pass_cut=false;

                  
                  ApplyCut(pass_cut,CAND1mTOFS[j]);
			
		}
            //cout << endl;
            
      }



}
