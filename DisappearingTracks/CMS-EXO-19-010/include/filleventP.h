#pragma once

#include <cassert>
#include <cstdlib>
#include <cstddef>
#include <type_traits>

#include <random>
#include <cfloat> // for DBL_MAX

//#include "time.h"
#include <chrono>

#include "heputil.h"
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

#include "Pythia8Plugins/CombineMatchingInput.h"
#include "Pythia8Plugins/aMCatNLOHooks.h"

#include "smearing.h"
#include <iomanip>

using namespace std;

HEP::P4 PVtoH(const Pythia8::Vec4 &pyv)
{
  HEP::P4 newv(pyv.px(), pyv.py(), pyv.pz(), pyv.e());
  return newv;
  //return HEP::P4(double pyv.px(), double pyv.py(), double pyv.pz(), double pyv.E());
}






bool isHSCP(Pythia8::Particle *particle, const Pythia8::Event &event)
//Return True (False) if the particle is (is not) a HSCP
{

	//Skip neutral particles
	if (!particle->isCharged()) return false;
	//Skip color charged particles or hadrons
	if (particle->colType() != 0) return false;
	if (particle->isHadron()) return false;
	//Skip light particles
	if (fabs(particle->m()) < 20.) return false;
	//Check if it is not an intermediate state (has itself as daughter)
	std::vector<int> daughters = particle->daughterListRecursive();
	for (int i = 0; i < daughters.size(); ++i){
		int idaughter = daughters[i];
		Pythia8::Particle daughterParticle = event[idaughter];
		if (daughterParticle.id() == particle->id()) return false;
	}

	return true;
}


bool decayBeforeEndHcal(Pythia8::Particle &particle)
//Return True (False) if the particle decays before (after) the hadronic calorimeter
{
	if (particle.isFinal()){return false;}
    Pythia8::Vec4 decayVertex = particle.vDec();
	if (decayVertex.pT() < 3.9e3 && abs(decayVertex.pz()) < 6.1e3){return true;}

	return false;
}

bool decayInsideAtlas(Pythia8::Particle &particle)
//Return True (False) if the particle decays inside the detector
{
	if (particle.isFinal()){return false;}
    Pythia8::Vec4 decayVertex = particle.vDec();
    if (decayVertex.pT() < 12e3 && abs(decayVertex.pz()) < 23e3){return true;}

	return false;
}

bool ischargedTrackHadron(Pythia8::Particle* particle, const Pythia8::Event &event)
{
  if(!particle->isHadron()) {return false;}
  if(!particle->isCharged()) { return false;}

  Pythia8::Vec4 prodVertex = particle->vProd();
  //if( (prodVertex.pT() > 1100.) || (abs(prodVertex.pz()) > 3100.)) {return false;}
  if( (prodVertex.pT() > 200.) || (abs(prodVertex.pz()) > 270.)) {return false;}
  
  if(particle->isFinal())
  {
    return true;
  }
  else
  {
        std::vector<int> daughters = particle->daughterListRecursive();
	      for (int i = 0; i < daughters.size(); i++){
		    int idaughter = daughters[i];
		      const Pythia8::Particle* daughterParticle = &event[idaughter];
		    if (daughterParticle->id() == particle->id()) return false;
	      }
        
        if( (particle->vDec().pT() < 300.0) ){return false;}
  }

  return true;
}




    inline bool fromHadron(int n, const Pythia8::Event &evt) {
      // Root particle is invalid
      if (n == 0) return false;
      const Pythia8::Particle p = evt[n];
      if (p.isHadron()) return true;
      if (p.isParton()) return false; // stop the walking at the end of the hadron level
      for (int m : p.motherList()) {
        if (fromHadron(m, evt)) return true;
      }
      return false;
    }


/*

bool registersinHCAL(Pythia8::Particle* prt)
{
  Pythia8::Vec4 pos = prt->vProd();
  

}


*/





void SortEvent(Pythia8::Event &event, HEP::Event &OutEvt, int npileups, std::vector<HEP::PileupEvent*> &pileups,double &Rjet, std::mt19937 &engine){

fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rjet);
std::vector <HEP::Jet*> outjets;
//std::vector <HEPUtils::Jet> outjets;
std::vector <fastjet::PseudoJet> inclusiveJets, jets;

std::vector <fastjet::PseudoJet> pileupJets, pJInputs; // for calculating rho
pJInputs.resize(0);



std::vector <fastjet::PseudoJet> fjInputs, bquark_momenta; //particles for applying jet clustering
    // Loop over particles in the event: kinematic distribution
fjInputs.resize(0);
bquark_momenta.resize(0);
//cout << "Starting search for jets" << endl;


//for(int mpart = 0; mpart< event.size(); mpart++)
for(int mpart = event.size()-1; mpart>=0 ; mpart--)
  {

    Pythia8::Particle* prt = &event[mpart];
    
    int pid=prt->id();
    int apid=abs(pid);
    
    
    if(prt->isFinal())
    {
     // This needs to either be added or deleted -- otherwise have a memory leak!
      HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
       

      if(prt->isHadron() && prt->isVisible())
         {
           if(prt->isCharged())
           {
             //OutEvt.add_charged_hadron(newpart->mom());
             Pythia8::Vec4 pos = prt->vProd();
             if (pos.pAbs() > 0.0001)
                  {
                    newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
                  }
             newpart->set_3Q(prt->chargeType());
             OutEvt.add_charged_hadron(newpart);

           }
           else
           {
             OutEvt.add_neutral_hadron(newpart->mom());
             delete newpart;
           }
         }
         else // not hadron
          {
            
              // look at the vertices
              	//if ( (*p)->production_vertex() ) { // i.e. does it have one: we want the position
                  
                  Pythia8::Vec4 pos = prt->vProd();


                  if (pos.pAbs() > 0.0001)
                  {
                    newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
                  }
                  

               // } 
            //newpart->set_prompt(); // sets it to not be prompt, otherwise the routine deletes it!!!
              if(apid == 11 || apid == 13 || apid == 22 ) // leptons or photons that *may* be "signal," or "background" and go into jets
              { 
                //if(!fromHadron(npart,event))
                //{
                  newpart->set_prompt(); // won't delete particle!

                  if(apid == 11 || apid ==13 )
                  {
                    if(pid > 0) {newpart->set_3Q(-3);} else {newpart->set_3Q(3);}

                  }


                  OutEvt.add_particle(newpart); 
                  
                  continue; // don't include in the jet
                /*}
                else
                {
                  if(apid == 22)
                  {
                    OutEvt.add_neutral_hadron(newpart->mom()); // photons contribute to calo energy like 'neutral hadrons'
                  }
                  else
                  {
                    OutEvt.add_charged_hadron(newpart->mom());     
                  }
                  
                }
                */
                
                 
              }  
              else
              {
                 // Stable particle, if it's a BSM particle let's find the charge
                 if(apid > 13)
                 {
                   newpart->set_3Q(prt->chargeType());
                 }                 
                 newpart->set_prompt(); 
               }
         


           OutEvt.add_particle(newpart); // particle is not prompt by default, so is deleted if it comes from a hadron
        }
      

      if( pid == 12 || pid == 14 || pid == 16 || pid == 18  ) // include electrons and muons, not neutrinos
      {
        continue;
      }
       if((pid > 1000000) && (pid < 9000000)) 
      {
        continue;
      }
      //HepMC::FourVector v4((*p)->momentum());
      Pythia8::Vec4 v4=prt->p();
      fastjet::PseudoJet particleTemp(v4.px(),v4.py(),v4.pz(),v4.e());
      fjInputs.push_back(particleTemp);
    }
    else
    {
      // cout << "is not final ";
       if(ischargedTrackHadron(prt,event))
       {
        
        Pythia8::Vec4 pos = prt->vProd();
        Pythia8::Vec4 vdec = prt->vDec();
        if(vdec.pT() < 33.0) continue; // check to make sure it doesn't decay too soon

        // Now must allocate or delete
        HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
        if (pos.pAbs() > 0.0001)
        {
          newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
        }
        newpart->set_decay(HEP::P4(vdec.px(),vdec.py(),vdec.pz(),vdec.e()));
        newpart->set_3Q(prt->chargeType());
    
        OutEvt.add_charged_hadron(newpart);
        continue;

       }


      /*
      if(apid == 15) // tau
      {
        bool islast = true;

        std::vector<int> daughters = prt.daughterList();
        for(int di = 0; di < daughters.size(); di ++)
        {
          if (event[daughters[di]].idAbs() == 15)
          {
            islast=false;
            break;
          }

        }
        
        if(islast)
        {
          HEP::Particle* newpart = new HEP::Particle(prt.px(),prt.py(),prt.pz(),prt.e(),pid);
          newpart->set_prompt();  
          Pythia8::Vec4 pos = prt.vProd();
          if (pos.pAbs() > 0.001)
          {
            newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
          }
          pos = prt.vDec();
          if (pos.pAbs() > 0.001)
          {
            newpart->set_decay(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
          }

          OutEvt.add_particle(newpart);
          continue;
        }
      }
    */
      
      if(isHSCP(prt,event))
      {
        //cout << "it is!" << prt->vDec().pT() << endl;
	//cout << prt->id() << " it is!" << prt->vDec().pAbs() << " : " << prt->tau() << ", " << prt->tau0() << endl;

        HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
        Pythia8::Vec4 pos = prt->vProd();

        if (pos.pAbs() > 0.001)
        {
          newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
        }
        else
        {
          newpart->set_prompt();  
        }
        //pos.clear();
        pos = prt->vDec();

        if (pos.pAbs() > 0.0001)
        {
          newpart->set_decay(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
        }
        newpart->set_3Q(prt->chargeType());
        OutEvt.add_HSCP(newpart);
        continue;
      }
      
      /*
     
      if(apid == 5)// use quarks as proxy, , (ought to use b hadrons, but whatever)
       {
    // ought to check to make sure they are part of the final process and not the incoming one, but whatever, I doubt it matters!
      // 
        //Make sure it's the last b quark
        bool islast = true;

        std::vector<int> daughters = prt.daughterList();
        for(int di = 0; di < daughters.size(); di ++)
        {
          if (event[daughters[di]].idAbs() == 5)
          {
            islast=false;
            break;
          }

        }

        if(islast)
        {
          Pythia8::Vec4 v4 = prt.p();
        
          fastjet::PseudoJet tpjet(v4.px(),v4.py(),v4.pz(),v4.e());

          bquark_momenta.push_back(tpjet);
        }
        }
        */
        //bquark_momenta.push_back(fastjet::PseudoJet(v4.px(),v4.py(),v4.pz(),v4.e()));
       
       //cout << "check done" << endl;
    }
  };

std::uniform_int_distribution<> rd_el(0,pileups.size()-1);

static std::normal_distribution<double> zdist(0,53.0);
static std::normal_distribution<double> tdist(0,160e-12);



fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
inclusiveJets = clustSeq.inclusive_jets(1.0); // argument is pTjet
jets    = sorted_by_pt(inclusiveJets);




pJInputs.clear();

/* Enable for pileup subtraction */
/*
for (auto const pJ : fjInputs)
{
  pJInputs.push_back(pJ);
}
*/

for(int jp=0; jp< npileups; jp++)
{
  int pileup_number=rd_el(engine);

  double dz=zdist(engine);
  double dt=tdist(engine);
  OutEvt.add_particles(pileups[pileup_number]->translate_particles(dz,dt));
  OutEvt.add_neutral_hadrons(pileups[pileup_number]->get_neutral_hadrons());
 
  OutEvt.add_charged_hadrons(pileups[pileup_number]->translate_charged_hadrons(dz,dt));


  /* Enable for pileup subtraction */
  /*
  // Add pileup to the computation of the median rho
  for(auto const chad : pileups[pileup_number]->get_charged_hadrons())
  {
    fastjet::PseudoJet particleTemp(chad->mom().px(),chad->mom().py(),chad->mom().pz(),chad->E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

  }

  for(auto const neut : pileups[pileup_number]->get_neutral_hadrons())
  {
    fastjet::PseudoJet particleTemp(neut->px(),neut->py(),neut->pz(),neut->E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

  }

  for(auto const vis : pileups[pileup_number]->get_particles())
  {
    fastjet::PseudoJet particleTemp(vis->mom().px(),vis->mom().py(),vis->mom().pz(),vis->mom().E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

  }
  */

 
}





//fastjet::ClusterSequence PJSeq(pJInputs, jetDef);
//pileupJets = PJSeq.inclusive_jets(1.0); // argument is pTjet
//GridMedianBackgroundEstimator bge(double max_rapidity,double requested_grid_spacing);
////bge.set_particles(event_particles);
//// the median of (pt/area) for grid cells, or for jets that pass the selection cut,
//// making use also of information on empty area in the event (in the jets case)
//rho = bge.rho();
//// an estimate of the fluctuations in theptdensity per unit√A,
//sigma = bge.sigma();

// 0.6 according to https://arxiv.org/pdf/1607.03663.pdf
// I'll use 0.5 since it's the cone size 

/*
fastjet::GridMedianBackgroundEstimator bge(4.7,0.5);
//bge.set_particles(pileupJets);
bge.set_particles(pJInputs);
OutEvt.Average_rho=bge.rho();
*/
pileupJets.clear();





// now do b tagging via method of looking if the delta_R of each jet to each initial b_quark is large enough
bool btag;
for(auto const& jjet: jets)
{
  btag=false;
  for( auto const& bquark : bquark_momenta)
  {
    if(bquark.delta_R(jjet) < Rjet) { btag = true; break;}
  }
  // now put into a ExHEPUtils jet with the appropriate tag

  // but apply the 77% efficiency in b-tagging:
  //smearbtag(&btag,jjet.perp(),engine);

  HEP::Jet* tjet = new HEP::Jet(jjet.px(),jjet.py(),jjet.pz(),jjet.E(),btag);

  outjets.push_back(tjet);
}



// now smear the jets once we have efficiencies for CMS
//smearJets(outjets,engine);
 


OutEvt.set_jets(outjets);

}          
          






Pythia8::Vec4 getMissingMomentum(Pythia8::Event &event)
{
  Pythia8::Vec4 missingETvec;
  std::vector<int> metaStableMothers;
    //First check for all final particles:
    for (int i = 0; i < event.size(); ++i) {
    	Pythia8::Particle ptc = event[i];
    	if (!ptc.isFinal()) continue;
    	//If particle has been produced inside or after the calorimeter, skip it and include its parent
      // These for ATLAS
    	//if (abs(ptc.vProd().pz()) > 6100. || ptc.vProd().pT() > 2300.){
        // These for CMS
        //if (abs(ptc.vProd().pz()) > 5500. || ptc.vProd().pT() > 2950.){
          // CMS muon chamber. Assuming here something invisible produced by something charged
          if (fabs(ptc.vProd().pz()) > 11000. || ptc.vProd().pT() > 7000.){
    		int mom = ptc.mother1();
    		//cout << "Adding mom " << mom << " from " << ptc.name() << " and index " << i << endl;
    		if(std::find(metaStableMothers.begin(), metaStableMothers.end(), mom) == metaStableMothers.end()) {
    			metaStableMothers.push_back(mom);
    		}
    		continue;
      	}
      
      
    	//Neutrinos are included in MET:
    	if (ptc.idAbs() == 12 || ptc.idAbs() == 14 || ptc.idAbs() == 16){
    		//cout << "Including: " << i << " " << ptc.name() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	if (ptc.isHadron()) continue;  //Ignore usual hadrons
    	if (ptc.idAbs() == 22) continue; //Ignore photons
      if (ptc.chargeType() !=0) continue; // Ignore anything charged!!
    	if (ptc.isVisible() && ptc.m() < 10.) continue; //Ignore light visible particles

        //Add to missing momentum
    	//cout << "Including: " << i << " " << ptc.name() << endl;
        missingETvec += ptc.p();
    }

    //Now check all meta-stable mothers:
    for (int imom = 0; imom < metaStableMothers.size(); ++imom) {
    	Pythia8::Particle ptc = event[metaStableMothers[imom]];
      
      
    	//Neutrinos are included in MET:
    	if (ptc.idAbs() == 12 || ptc.idAbs() == 14 || ptc.idAbs() == 16){
    		//cout << "Including mom: " << metaStableMothers[imom] << " " << ptc.id() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	if (ptc.isHadron()) continue;  //Ignore usual hadrons
    	if (ptc.idAbs() == 22) continue; //Ignore photons
       if (ptc.chargeType() !=0) continue; // Ignore anything charged!!
    	if (ptc.isVisible() && ptc.m() < 10.) continue; //Ignore light visible particles

        //Add to missing momentum
    	//cout << "Including mom: " << metaStableMothers[imom] << " " << ptc.id() << endl;
        missingETvec += ptc.p();
    }

    return missingETvec;


}

/*
void calcMissingMomentum(HEP::Event &OutEvt, std::mt19937 &engine)
{
  // look at the invisibles. Bit stupid I know but whatever. 
  HEP::P4 pmiss;
  pmiss.clear();
  for(HEP::Particle* invisible : OutEvt.invisible_particles())
  {
    pmiss+=invisible->mom();
  }

  OutEvt.set_missingmom(pmiss);

  // smear here ...
}
*/

HEP::PileupEvent* sortPileupEvent(Pythia8::Event &event)
{
  HEP::PileupEvent* newPuE = new HEP::PileupEvent();


  for(int npart = 0; npart< event.size(); npart++)
  {
    Pythia8::Particle* prt = &event[npart];
    int pid=prt->id();
    int apid=abs(pid);
    // check if is final and produced inside tracker or if it decays outside enough (not yet)
    //cout << pid << endl;
    if(prt->isFinal())
    {
      HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
      // All particles other than invisibles and muons are jet constituents. Bollocks, really we only want photons
      //if( pid == 12 || pid == 14 || pid == 16 || pid == 18 || pid == 13 || pid == 11)
      

      if(prt->isHadron() && prt->isVisible())
         {
           if(prt->isCharged())
           {
             //OutEvt.add_charged_hadron(newpart->mom());
             Pythia8::Vec4 pos = prt->vProd();
             if (pos.pAbs() > 0.0001)
                  {
                    newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
                  }
             newpart->set_3Q(prt->chargeType());
             
             //newPuE->add_hadron_vertex(newpart);
             if(newpart->pT() > 1.0)
             {
              newPuE->add_charged_hadron(newpart);
             }
             else
             {
               delete newpart;
             }
           }
           else
           {
             newPuE->add_neutral_hadron(newpart->mom());
             delete newpart;
           }
         }
         else // not hadron
          {
            
              // look at the vertices
              	//if ( (*p)->production_vertex() ) { // i.e. does it have one: we want the position
                  //HepMC::FourVector pos =  (*p)->production_vertex()->position();
                  Pythia8::Vec4 pos = prt->vProd();


                  if (pos.pAbs() > 0.0001)
                  {
                    newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
                  }
                  newpart->set_prompt(); 
                  newpart->set_3Q(prt->chargeType());
                  newPuE->add_particle(newpart);
               // } 
            
          }   
      // don't bother with jet constituents for this analysis
      /*
      if( pid == 12 || pid == 14 || pid == 16 || pid == 18  ) // include electrons and muons, not neutrinos
      {
        continue;
      }
       if((pid > 1000000) && (pid < 9000000)) 
      {
        continue;
      }
      //HepMC::FourVector v4((*p)->momentum());
      Pythia8::Vec4 v4=prt.p();
      fastjet::PseudoJet particleTemp(v4.px(),v4.py(),v4.pz(),v4.e());
      fjInputs.push_back(particleTemp);
      */
    }
    else
    {
      // cout << "is not final ";
       if(ischargedTrackHadron(prt,event))
       {
        
        Pythia8::Vec4 pos = prt->vProd();
        Pythia8::Vec4 vdec = prt->vDec();
        if(vdec.pT() < 33.0) continue; // check to make sure it doesn't decay too soon
        HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
              if (pos.pAbs() > 0.0001)
                  {
                    newpart->set_prod(HEP::P4(pos.px(),pos.py(),pos.pz(),pos.e()));
                  }
              newpart->set_decay(HEP::P4(vdec.px(),vdec.py(),vdec.pz(),vdec.e()));
             newpart->set_3Q(prt->chargeType());
            // cout << " adding charged hadron ";
             //newPuE->add_hadron_vertex(newpart);
             newPuE->add_charged_hadron(newpart);


       }
    }


  }

  return newPuE;

}




void filleventPYTHIA(HEP::Event &OutEvt, Pythia8::Event &evt, double weight, int npileups, std::vector<HEP::PileupEvent*> &pileups, std::mt19937 &engine)
{
// Function for filling a heputils event from the HepMC one, doing the smearing etc.

// First let's set the weight

OutEvt.set_weight(weight);

double dR=0.4;
std::vector<HEP::Jet*> tjets;

SortEvent(evt,OutEvt,npileups,pileups,dR,engine);
// fill up the jets
//OutEvt.set_jets(tjets);

// fill up final state particles, compute MET



OutEvt.set_missingmom(PVtoH(getMissingMomentum(evt)));

//OutEvt.set_missingmom(smearMET(OutEvt.missingmom(),engine));

// add smearing once we have efficiencies for CMS
/*
filterEfficiency(OutEvt.jets(),JetEff,engine);


filterEfficiency(OutEvt.electrons(),ElectronEff,engine);
filterEfficiency(OutEvt.muons(),MuonEff,engine);

smearElectronEnergy(OutEvt.electrons(),engine);
smearMuonMomentum(OutEvt.muons(),engine);
*/



}
