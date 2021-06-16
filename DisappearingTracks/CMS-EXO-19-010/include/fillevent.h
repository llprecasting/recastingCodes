#pragma once

#include <cassert>
#include <cstdlib>
#include <cstddef>
#include <type_traits>

#include <random>
#include <cfloat> // for DBL_MAX


#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"


#include "heputil.h"

#include "smearing.h"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

using namespace std;

bool isfinal(const HepMC::GenParticle* p )
{
//if ( (!(p->end_vertex())) && (p->status()==1) ) return true;
//if(p->status() ==3) return false;

if ( (!(p->end_vertex()))) 
{
  if(p->status() != 1) 
  {
    cout << "Particle with no end vertex and status" << p->status() << endl;
    return false;
  }
  return true;
}

return false;

};


/*

HepMC: perp() gives pT, rho() gives 3-vector length



*/


bool isDT(const HepMC::GenParticle* p, const HepMC::GenEvent *event)
//Return True (False) if the particle is (is not) a HSCP
{
int pid = p->pdg_id();
 int pch = PDG::get3Q(pid);

	//Skip neutral particles
	if (pch == 0) return false;
	
  
  // Care if it's a hadron! The metastable charged hadrons will be included elsewhere ...
  if(PDG::isHadron(pid)) return false;
	
  
  //Skip light particles
	//if (fabs(particle->m()) < 20.) return false;
  if(p->momentum().m() < 20.) return false;

  

	//Check if it is not an intermediate state (has itself as daughter)
  for ( HepMC::GenVertex::particle_iterator des = p->end_vertex()->particles_begin(HepMC::descendants); des != p->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if ((*des)->pdg_id() == pid)
          {
            return false;
 
          }
		    }

  
  // check produced and decays inside
  if(!p->end_vertex()) return false;
  if(!p->production_vertex()) return false;

  HepMC::FourVector pos =  p->production_vertex()->position();
  HepMC::FourVector dec = p->end_vertex()->position();

  double pperp = pos.perp();
  double dperp = dec.perp();
  double pz = pos.pz();
  double dz = dec.pz();
  
  // check if track is long enough
  if(((dperp - pperp) < 3.0 ) && (fabs(dz-pz) < 3.0)) return false;
  
  // check whether it shows up as a track
  //if((pperp > 1100.0) || (dperp < 160.0) || (fabs(pz) > 2700.0) || (fabs(dz) > 2700.0)) return false;
  if((pperp > 1100.0) || (dperp < 160.0) || (fabs(pz) > 2700.0)) return false;
  //cout << "so far so good " << pid <<  " " << pz << ", " << dz << "," << dperp << std::endl;
	//std::vector<int> daughters = particle->daughterListRecursive();
	//for (int i = 0; i < daughters.size(); ++i){
	//	int idaughter = daughters[i];
	//	Pythia8::Particle daughterParticle = event[idaughter];
	//	if (daughterParticle.id() == particle->id()) return false;
	//}

	return true;
}

bool isMetaStableHadron(const HepMC::GenParticle* p, const HepMC::GenEvent *event)
//Return True (False) if the particle is (is not) a HSCP
{
int pid = p->pdg_id();
// int pch = PDG::get3Q(pid);
 //Skip neutral particles
//	if (pch == 0) return false;
 // skip hadrons
if(!(PDG::isHadron(pid))) return false;
	
	

  
	//Skip light particles
	//if (fabs(particle->m()) < 20.) return false;
  if(p->momentum().m() < 20.) return false;
  
  // We're interested in the tracker
  //if( (p->production_vertex()->position().perp() > 1100.) || (fabs(p->production_vertex()->position().pz()) > 3100.)) {return false;}
  // Let's say must be produced in the pixel detector, decays somewhere in the tracker or outside. That way, we include it in
  // the isolation calculation only once

  double prodT=p->production_vertex()->position().perp();
  double prodZ=p->production_vertex()->position().pz();

  if( (prodT > 200.) || (fabs(prodZ) > 270.)) {return false;}

  


	//Check if it is not an intermediate state (has itself as daughter)
  for ( HepMC::GenVertex::particle_iterator des = p->end_vertex()->particles_begin(HepMC::descendants); des != p->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if ((*des)->pdg_id() == pid)
          {
            return false;
 
          }
		    }
  double decT=p->end_vertex()->position().perp();
  double decZ=p->end_vertex()->position().pz();
  
  //if( (p->end_vertex()->position().perp() < 33.0)){return false;}
  if( (decT < 300.0) ){return false;}

  

	//std::vector<int> daughters = particle->daughterListRecursive();
	//for (int i = 0; i < daughters.size(); ++i){
	//	int idaughter = daughters[i];
	//	Pythia8::Particle daughterParticle = event[idaughter];
	//	if (daughterParticle.id() == particle->id()) return false;
	//}

	return true;
}





bool isHSCP(const HepMC::GenParticle* p, const HepMC::GenEvent *event)
//Return True (False) if the particle is (is not) a HSCP
{
int pid = p->pdg_id();
 int pch = PDG::get3Q(pid);
 //Skip neutral particles
	if (pch == 0) return false;
 // skip hadrons
if(PDG::isHadron(pid)) return false;
	
	

  
	//Skip light particles
	//if (fabs(particle->m()) < 20.) return false;
  if(p->momentum().m() < 20.) return false;

  //if( (prodVertex.pT() > 1100.) || (abs(prodVertex.pz()) > 3100.)) {return false;}
 
	//Check if it is not an intermediate state (has itself as daughter)
  for ( HepMC::GenVertex::particle_iterator des = p->end_vertex()->particles_begin(HepMC::descendants); des != p->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if ((*des)->pdg_id() == pid)
          {
            return false;
 
          }
		    }

	//std::vector<int> daughters = particle->daughterListRecursive();
	//for (int i = 0; i < daughters.size(); ++i){
	//	int idaughter = daughters[i];
	//	Pythia8::Particle daughterParticle = event[idaughter];
	//	if (daughterParticle.id() == particle->id()) return false;
	//}

	return true;
}

/*
inline bool fromHadron(int n, const Pythia8::Event& evt) {
      // Root particle is invalid
      if (n == 0) return false;
      const Pythia8::Particle& p = evt[n];
      if (p.isHadron()) return true;
      if (p.isParton()) return false; // stop the walking at the end of the hadron level
      for (int m : p.motherList()) {
        if (fromHadron(m, evt)) return true;
      }
      return false;
    }

bool isfinal(const HepMC::GenParticle* p )
{
if ( !p->end_vertex() && p->status()==1 ) return true;
return false;

};
*/

/*
double d0Calc(const HepMC::FourVector xv, const HepMC::FourVector p)
//Compute the transverse impact parameter for a particle with momentum p created at vertex xv
//assuming no magnetic field
// Look at x-y plane and find clsost distance: this is |xv| sin varphi, where the angle is between the vertex
// and the momentum
// this is just given by xv cross phat, where each is represented as a 2d vector
{
    double calpha = p.px()/p.perp();
    double salpha = p.py()/p.perp();

    double a0 = xv.px()*salpha - xv.py()*calpha;

    return fabs(a0);
}

double dzCalc(const HepMC::FourVector xv, const HepMC::FourVector p)
//Compute the longitudinal impact parameter for a particle with momentum p created at vertex xv
//assuming no magnetic field
// dz is the z value at the point of closest d0.
// v_min = xv - lambda \hat{p}
// lambda = |xv| cos varphi = xv dot \hat{p} (in 2d)
{
    //double calpha = p.px()/p.pT();
    //double salpha = p.py()/p.pT();

    //double a0 = xv.px()*salpha - xv.py()*calpha;
    double lambda= (p.px()*xv.px()+p.py()*xv.py())/p.perp();
    double zmin = xv.pz()-lambda*p.pz()/p.perp();
    return fabs(zmin);
}

double dz_sintheta(const HepMC::FourVector xv, const HepMC::FourVector p)
//Compute the longitudinal impact parameter for a particle with momentum p created at vertex xv
{
    //double calpha = p.px()/p.pT();
    //double salpha = p.py()/p.pT();

    //double a0 = xv.px()*salpha - xv.py()*calpha;
    if(p.rho() == 0.0) return 0.0;
    double dz = dzCalc(xv,p);

    return fabs(dz*p.perp()/p.rho());
}
*/

void DebugMyElectron(HepMC::GenEvent *event,const HepMC::GenParticle* p )
{

const HepMC::GenParticle* q = new HepMC::GenParticle;
q=p;
cout << "Diagnosing" << endl;
int pid = p->pdg_id();
bool foundGenuineMother=false;
while((!foundGenuineMother) && (q->production_vertex()))
{
  for ( HepMC::GenVertex::particle_iterator mother= q->production_vertex()->particles_begin(HepMC::parents);
			  mother != q->production_vertex()->particles_end(HepMC::parents); ++mother)
  {
    cout << "mother id: " << (*mother)->pdg_id() << ", ";
    if((*mother)->pdg_id() == pid)
    {
      q=*mother;
      break;
    }
    else
    {
      foundGenuineMother=true;
      break;
    }
  }
}
cout << endl;
}


bool isFromHadron(HepMC::GenEvent *event,const HepMC::GenParticle* p )
{
  // follow the chain up until we get to a vertex where there are no mothers equal to p
  //if ( !(*p)->production_vertex() ) return false; // weird case here, the particle has no production vertex???
  
  bool foundGenuineMother=false;
  int pid = p->pdg_id();
  const HepMC::GenParticle* q = new HepMC::GenParticle;

  int mpid;
  bool ismotherhadron;
  q = p;
  while((!foundGenuineMother) && (q->production_vertex()))
  { 
  foundGenuineMother= true;
  for ( HepMC::GenVertex::particle_iterator mother= q->production_vertex()->particles_begin(HepMC::parents);
			  mother != q->production_vertex()->particles_end(HepMC::parents); 
			  ++mother )
        {
          if((*mother)->pdg_id() == pid)
          {
            q = *mother; // the particle iterator points to the pointer to the particle
            foundGenuineMother=false;
            break;
          }
          else
          {
            mpid=(*mother)->pdg_id();
            ismotherhadron=PDG::isHadron(mpid);
            if(ismotherhadron)
            {
              foundGenuineMother=true;
              q = *mother;
              break;
            }
          }
          
        }
  };

  // here we've found the genuine mother
  if(ismotherhadron)
  {
    cout << " Hadronic mother " << mpid << endl;
    return true;
  }

  /*
  else
  {
    cout << " non-hadronic mother(s): " ;
    for ( HepMC::GenVertex::particle_iterator mother= q->production_vertex()->particles_begin(HepMC::parents);
			  mother != q->production_vertex()->particles_end(HepMC::parents); 
			  ++mother )
        {
          cout << (*mother)->pdg_id() << ", ";
        }
    
    cout << endl;
  }
  */

  return false; //ismotherhadron;

} 




void SortEvent(HepMC::GenEvent *event, HEP::Event &OutEvt, int npileups, std::vector<HEP::PileupEvent*> &pileups,double &Rjet, std::mt19937 &engine){

fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rjet);
std::vector <HEP::Jet*> outjets;
//std::vector <HEPUtils::Jet> outjets;
std::vector <fastjet::PseudoJet> inclusiveJets, jets;
std::vector <fastjet::PseudoJet> pileupJets, pJInputs; // for calculating rho

std::vector <fastjet::PseudoJet> fjInputs, bquark_momenta; //particles for applying jet clustering
    // Loop over particles in the event: kinematic distribution
fjInputs.resize(0);
pJInputs.resize(0);
bquark_momenta.resize(0);
//cout << "Starting search for jets" << endl;
HepMC::GenEvent::particle_iterator p = event->particles_begin();

//HEP::P4 p4miss(0.0,0.0,0.0,0.0);

while( p != event->particles_end())
  {
   
    int pid=(*p)->pdg_id();
    int apid=abs(pid);
    
    //cout << pid << endl;
    if( isfinal(*p))
    {
      HEP::Particle* newpart = new HEP::Particle((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e(),pid);
      // All particles other than invisibles and muons are jet constituents. Bollecks, really we only want photons
      //if( pid == 12 || pid == 14 || pid == 16 || pid == 18 || pid == 13 || pid == 11)
      if ( (*p)->production_vertex() ) 
        {
        
              HepMC::FourVector pos =  (*p)->production_vertex()->position();
                  if (pos.rho() > 0.0001)
                  {
                    // check if it's actually produced inside the detector
                    if((pos.perp() > 12e3) || (abs(pos.z()) > 23e3 )) { p++; delete newpart; continue;}
                    newpart->set_prod(HEP::P4(pos.x(),pos.y(),pos.z(),pos.t()));
                  }
        }

      if(PDG::isHadron(pid))
         {
           //p4miss-=newpart->mom();
           int pcharge=PDG::get3Q(pid);
           //if(PDG::isCharged(pid))
           if(pcharge != 0)
           {
             //OutEvt.add_charged_hadron((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());

             
              newpart->set_3Q(pcharge);
              OutEvt.add_charged_hadron(newpart);
           }
           else
           {
             OutEvt.add_neutral_hadron((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
             delete newpart;
           }
         }
         else
          {
            newpart->set_prompt(); // sets it to be prompt, otherwise the routine deletes it!!!

            if(apid == 22)
              {
                newpart->set_3Q(0);
               
              }

            if((apid == 11) || (apid == 13))
            { 
              // Some debug code when I was looking at displaced vertex decays 
              // Turned out you can get what looks like a four-electron vertex when an emitted photon
              // produced an electron-positron pair in the same place as the main decay
              /*
               if((newpart->prod_vertex().pT()) > 1530 && (newpart->prod_vertex().pT() < 1532))
               {
                 cout << "Dodgy electron " << newpart->mom() << endl;
                 DebugMyElectron(event,*p);
               }
               */
                if(pid > 0)
                {
                  newpart->set_3Q(-3);

                }
                else
                {
                  newpart->set_3Q(3);
                }
            }
            else
            {
              if(apid > 13) // care about heavy stable charged particles 
              {
                newpart->set_3Q(PDG::get3Q(pid));
              }
            }      

            OutEvt.add_particle(newpart);
          }
       

      if( pid == 12 || pid == 14 || pid == 16 || pid == 18  ) // include electrons and muons in jets
      {
        
        
        p++;
        continue;
      }
      /*
      if(PDG::get3Q(pid) != 0)
        {
          p4miss-=newpart->mom();
        }
        */
       if((pid > 1000000) && (pid < 9000000)) 
      {
        
        p++;
        continue;
      }
      HepMC::FourVector v4((*p)->momentum());
      //cout << "v4x: " << v4.px() << endl;


      fastjet::PseudoJet particleTemp(v4.px(),v4.py(),v4.pz(),v4.e());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      fjInputs.push_back(particleTemp);
      
        p++;
      //constituent_pids.push_back(apid);
      //if(njetparts != constituent_pids.size())
      //{
      //  cout << "Error in size of pids in jets!" << endl;
      //}
    }
    else
    {
      //cout << "check here" << endl;
      if(!((*p)->end_vertex())) {  p++; continue;}
      
      if(isDT(*p,event)) 
      {
  //      cout << "DT: " << pid  << std::endl;
        HEP::Particle* newpart = new HEP::Particle((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e(),pid);
        HepMC::FourVector pos =  (*p)->production_vertex()->position();           
        newpart->set_prod(HEP::P4(pos.x(),pos.y(),pos.z(),pos.t()));
        HepMC::FourVector dec =  (*p)->end_vertex()->position();
        newpart->set_decay(HEP::P4(dec.x(),dec.y(),dec.z(),dec.t()));
        newpart->set_3Q(PDG::get3Q(pid));
        newpart->set_prompt();
        OutEvt.add_particle(newpart);  
        p++; continue;        
      }

      if(isMetaStableHadron(*p,event))
      {
        HEP::Particle* newpart = new HEP::Particle((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e(),pid);
        HepMC::FourVector pos =  (*p)->production_vertex()->position();           
        newpart->set_prod(HEP::P4(pos.x(),pos.y(),pos.z(),pos.t()));
        HepMC::FourVector dec =  (*p)->end_vertex()->position();
       
        newpart->set_decay(HEP::P4(dec.x(),dec.y(),dec.z(),dec.t()));
        int ch = PDG::get3Q(pid);
        newpart->set_3Q(ch);
        newpart->set_prompt();  
        if(ch !=0)
        {
          
          OutEvt.add_charged_hadron(newpart);
          //OutEvt.add_HSCP(newpart);

        }
        else
        {
          delete newpart;
        }
        
      }

      // don't care about b-tagging
      /*
      if(((*p)->end_vertex()) && (apid == 5))// use quarks as proxy, , (ought to use b hadrons, but whatever)
       {
    // ought to check to make sure they are part of the final process and not the incoming one, but whatever, I doubt it matters!
      // 
        //Make sure it's the last b quark
        bool islast = true;
        for ( HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()->particles_begin(HepMC::descendants); des != (*p)->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if (abs((*des)->pdg_id()) == 5)
          {
            islast=false;
            break;
          }
		    }
        if(islast)
        {
          HepMC::FourVector v4((*p)->momentum());
        
          fastjet::PseudoJet tpjet(v4.px(),v4.py(),v4.pz(),v4.e());

          bquark_momenta.push_back(tpjet);
        }
        
        //bquark_momenta.push_back(fastjet::PseudoJet(v4.px(),v4.py(),v4.pz(),v4.e()));
       }  
       */

       
      p++;

       //cout << "check done" << endl;
    }
  
  };

/*std::cout << "Found " << bquark_momenta.size() << " b quarks" << endl;
for( auto bm : bquark_momenta)
{
  cout << bm.perp() << ", ";

}
cout << endl;
*/

/*
for (int i = 0; i < event.size(); ++i) {
        // Require visible particles inside detector.
if (! isfinal(event[i]) || event[i].isLepton() || abs(event[i].id()) == 1000024 || abs(event[i].id()) == 1000022
      || event[i].id() == 23
      || abs(event[i].id()) == 24){continue;}
        fastjet::PseudoJet particleTemp(event[i].px(),event[i].py(),event[i].pz(),event[i].e());
        particleTemp.set_user_index(i);
        fjInputs.push_back(particleTemp);
    }

*/    
fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
inclusiveJets = clustSeq.inclusive_jets(1.0); // argument is pTjet
    jets    = sorted_by_pt(inclusiveJets);

std::uniform_int_distribution<> rd_el(0,pileups.size()-1);

static std::normal_distribution<double> zdist(0,53.0);
static std::normal_distribution<double> tdist(0,160e-12);

//auto starttime = std::chrono::high_resolution_clock::now();
//
//time_t starttime = time(NULL);

for(int jp=0; jp< npileups; jp++)
{
  int pileup_number=rd_el(engine);
  //cout << "selecting pileup number " << pileup_number << std::endl;
  //for(auto nhad : pileups[pileup_number]->get_neutral_hadrons())
 // {
  //  OutEvt.add_neutral_hadron(nhad)
 // }
  double dz=zdist(engine);
  double dt=tdist(engine);
  OutEvt.add_particles(pileups[pileup_number]->translate_particles(dz,dt));
  OutEvt.add_neutral_hadrons(pileups[pileup_number]->get_neutral_hadrons());
  //OutEvt.add_charged_hadrons(pileups[pileup_number]->get_charged_hadrons());
  OutEvt.add_charged_hadrons(pileups[pileup_number]->translate_charged_hadrons(dz,dt));
 // OutEvt.add_hadron_vertices(pileups[pileup_number]->translate_hadron_vertices(dz,dt));
}

// now to compute Jets including pileup
pJInputs.clear();
for(auto const chad : OutEvt.get_charged_hadrons())
{
  fastjet::PseudoJet particleTemp(chad->mom().px(),chad->mom().py(),chad->mom().pz(),chad->E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

}

for(auto const neut : OutEvt.get_neutral_hadrons())
{
  fastjet::PseudoJet particleTemp(neut->px(),neut->py(),neut->pz(),neut->E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

}

for(auto const vis : OutEvt.visible_particles())
{
  fastjet::PseudoJet particleTemp(vis->mom().px(),vis->mom().py(),vis->mom().pz(),vis->mom().E());
      //particleTemp.set_user_index(njetparts);
      //njetparts++;
      pJInputs.push_back(particleTemp);

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
fastjet::GridMedianBackgroundEstimator bge(4.7,0.5);
//bge.set_particles(pileupJets);
bge.set_particles(pJInputs);
OutEvt.Average_rho=bge.rho();

pileupJets.clear();

//auto endtime = std::chrono::high_resolution_clock::now();
//std::chrono::duration<double> timediff = std::chrono::duration_cast<std::chrono::duration<double>>(endtime-starttime);
//double timediff = difftime(endtime, starttime) * 1000.0;

/*
cout.precision(10);
cout << "Filling pileup took " << timediff.count() << std::endl;
cout << "Event now has " << OutEvt.get_charged_hadrons().size() << " charged hadrons, "  << OutEvt.get_neutral_hadrons().size() << " neutral hadrons, "  << OutEvt.visible_particles().size() << " particles " << std::endl;
cout << "Now filling event of n particles: " << event.size() << std::endl ;
*/




// now do b tagging via method of looking if the delta_R of each jet to each initial b_quark is large enough
//int btagged=0;
//std::uniform_real_distribution<double> rd(0.0,1.0);
/*
for (auto n : electron_numbers)
{
  cout << n << ", " ;
}
cout << endl;
*/

//cout << "Found " << bquark_momenta.size() << "bquarks" << endl;
//cout << "Found " << candidate_electrons.size() << " electron candidates" << endl;
//bool btag;




//cout << "smearing jets" << endl;
// now smear the jets
//smearJets(outjets,engine);
 
for(auto jet : jets)
{
  HEP::Jet* tjet = new HEP::Jet(jet.px(),jet.py(),jet.pz(),jet.E());
  outjets.push_back(tjet);

}

OutEvt.set_jets(outjets);
// now return them

//return outjets;


}


void calcMissingMomentum(HEP::Event &OutEvt, std::mt19937 &engine)
{
  // look at the invisibles. Bit stupid I know but whatever. 
  HEP::P4 pmiss;
  pmiss.clear();

  /*
  if (abs(ptc.vProd().pz()) > 6100. || ptc.vProd().pT() > 2300.)

  */
    
  for(HEP::Particle* invisible : OutEvt.invisible_particles())
  {
   // std::cout << fabs(invisible->prod_vertex().pz()) << ", " << invisible->prod_vertex().pT() << std::endl;
    //if((fabs(invisible->prod_vertex().pz()) > 6100.) || (invisible->prod_vertex().pT() > 2300.)) continue;
    //if((fabs(invisible->prod_vertex().pz()) > 3150.) || (invisible->prod_vertex().pT() > 1200.)) continue;
    pmiss+=invisible->mom();
  }
  
  /*
  for(HEP::Particle* invisible : OutEvt.HSCPs())
  {
    //if((fabs(invisible->prod_vertex().pz()) > 6100.) || (invisible->prod_vertex().pT() > 2300.)) continue;
   // if((fabs(invisible->decay_vertex().pz()) < 6100.) || (invisible->decay_vertex().pT() < 2300.)) continue;
    pmiss+=invisible->mom();
  }
  */
  OutEvt.set_missingmom(pmiss);

  // smear here ...
}



void anotherCalcMissingMomentum(HepMC::GenEvent *event,HEP::Event &OutEvt)
{

  HEP::P4 pmiss;
  pmiss.clear();

  HepMC::GenEvent::particle_iterator p = event->particles_begin();
  //HepMC::FourVector HMCpmiss;
//HEP::P4 p4miss(0.0,0.0,0.0,0.0);

while( p != event->particles_end())
  {
   
    int pid=(*p)->pdg_id();
    int apid=abs(pid);
    //HEP::P4 vprod(0.0,0.0,0.0,0.0);
    double prodZ=0.0;
    double prodT=0.0;
    if ( (*p)->production_vertex() ) 
        {
          HepMC::FourVector pos =  (*p)->production_vertex()->position();
          if (pos.rho() > 0.0001)
                  {
                    prodT=pos.perp();
                    prodZ=fabs(pos.z());
                   
                    // If produced outside the detector then it is no use to us
                     if((prodT > 7e3) || (prodZ > 11e3 )) { p++;  continue;}
                    
                    

                  }

        }

    if( isfinal(*p))
    {
      // Need it to be produced before end of HCAL, unless it is produced
      // with something like a muon track (ignore that case)
      //if((prodT > 2.95e3) || (prodZ > 5.5e3 )) { p++;  continue;}
      if(PDG::isInvisible(pid))
      {
        pmiss+=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());

      }




    }
    else
    {
      // want the last *visible* particle that decays to an invisible outside the HCAL (but still inside the detector)
      // NO! Actually we want the last *invisible* metastable particle that decays outside
      if(!PDG::isInvisible(pid)) {p++; continue;} 
      // see if it decays to an invisible particle beyond the calo
       
       double decT = (*p)->end_vertex()->position().perp();
       double decz = fabs((*p)->end_vertex()->position().pz());
        // Decay must be outside the detector, we already know it's produced inside

        // If it decays outside the muon system then it doesn't contribute to pmiss!!
       //if((decT > 7e3 ) || (decz > 11e3)) { p++; continue;}

        // If it decays in the HCAL then it doesn't contribute, but we've already checked it's in the detector
        //if((decT < 2.95e3 ) && (decz < 5.5e3)) { p++; continue;}
	// change to the muon system; if it decays inside then it is no good, it has to match above
	// Here we are only interested in invisibles produced in the detector that escape
	if((decT < 7.0e3 ) && (decz < 11.0e3)) { p++; continue;}



       // 
       bool islast = true;
       //bool invisibleproduct = false;
	
      for ( HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()->particles_begin(HepMC::descendants); des != (*p)->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if ((*des)->pdg_id() == pid)
          {
            islast=false;
            break;
          }
	 /*	
          else
          {
            if(PDG::isInvisible((*des)->pdg_id())) 
            {
              invisibleproduct=true;
              
            }
          }
	  */
        }

        if(!islast) {p++; continue;}
        //if(invisibleproduct)
        //{
          pmiss+=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
        //}

    }

    p++;



  }
  OutEvt.set_missingmom(pmiss);
}


void betterCalcMissingMomentum(HepMC::GenEvent *event,HEP::Event &OutEvt)
{
  HEP::P4 pmiss;
  pmiss.clear();

  HepMC::GenEvent::particle_iterator p = event->particles_begin();
  //HepMC::FourVector HMCpmiss;
//HEP::P4 p4miss(0.0,0.0,0.0,0.0);

while( p != event->particles_end())
  {
   
    int pid=(*p)->pdg_id();
    int apid=abs(pid);
    //HEP::P4 vprod(0.0,0.0,0.0,0.0);
    double prodZ=0.0;
    double prodT=0.0;
    if ( (*p)->production_vertex() ) 
        {
          HepMC::FourVector pos =  (*p)->production_vertex()->position();
          if (pos.rho() > 0.0001)
                  {
                    prodT=pos.perp();
                    prodZ=fabs(pos.z());
                    if((prodT > 12e3) || (prodZ > 23e3 )) { p++;  continue;}

                    

                  }

        }


    
    //cout << pid << endl;
    if( isfinal(*p))
    {

      if(PDG::isHadron(pid))
         {
           // make sure produced before end of HCAL
           if((prodT > 3000.0)|| (prodZ > 5000.0)) {p++; continue;}
          pmiss-=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
         }
         else
         {
            if( (apid == 22) || (PDG::get3Q(pid) !=0 ) )
            {
              if((prodT > 1700.0)|| (prodZ > 4000.0)) {p++; continue;}
              pmiss-=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
            }
         }

    }
    else
    {
      // check last of its kind
      if((!(*p)->end_vertex()) || (!(*p)->production_vertex())) {p++; continue;}

      bool islast = true;
      for ( HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()->particles_begin(HepMC::descendants); des != (*p)->end_vertex()->particles_end(HepMC::descendants);
			  ++des ) {
          if ((*des)->pdg_id() == pid)
          {
            islast=false;
            break;
          }
        }

        if(!islast) { p++; continue;}

      double decZ= fabs((*p)->end_vertex()->position().pz());
      double decT=(*p)->end_vertex()->position().perp();

      if(PDG::isHadron(pid))
      {
        if((prodT > 3000.0)|| (prodZ > 5000.0) || (decZ < 4e3 ) ||(decT < 1.77e3)) {p++; continue;}
        pmiss-=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
      }
      else
      {
        //if( (apid == 22) || (apid == 11) )
        if(PDG::isCharged(pid))
            {
              // treat as a muon
              //if((prodT > 1700.0)|| (prodZ > 4000.0) || (decZ < 3.15e3 ) ||(decT < 1.2e3)) {p++; continue;}
              if( (prodT > 1.77e3) || (prodZ > 4000.0) || (decZ < 1.2e4) || (decT < 7.0e3)){ p++; continue;}
              pmiss-=HEP::P4((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
            }
      }

    }

    
  p++;
  }
  
  //HEP::P4 pmiss(HMCpmiss.px(),HMCpmiss.py(),HMCpmiss.pz(),HMCpmiss.E());
  OutEvt.set_missingmom(pmiss);


}





void fillevent(HEP::Event &OutEvt, HepMC::GenEvent* evt, int npileups, std::vector<HEP::PileupEvent*> &pileups, std::mt19937 &engine)
{
// Function for filling a heputils event from the HepMC one, doing the smearing etc.

// First let's set the weight
// this needs correcting
//OutEvt.set_weights(evt->weights());

OutEvt.set_weight(evt->weights()[0]);
//HEPUtils::Event OutEvt;
double dR=0.4;
//std::vector<HEP::Jet*> tjets;

//cout << "Getting jets " << endl;

//cout << "sorting event" << endl;
SortEvent(evt,OutEvt,npileups,pileups,dR,engine);
// fill up the jets
//OutEvt.set_jets(tjets);

// fill up final state particles, compute MET

//cout << "filling leptons" << endl;

//betterCalcMissingMomentum(evt, OutEvt);
anotherCalcMissingMomentum(evt, OutEvt);
//calcMissingMomentum(OutEvt, engine);

//cout << "Before returning we have " << OutEvt.electrons().size() << " electrons and " << OutEvt.muons().size() << " muons" << endl;
//cout << "Before smearing we have " << OutEvt.electrons().size() << " electrons and " << OutEvt.muons().size() << " muons" << endl;


/*
filterEfficiency(OutEvt.jets(),JetEff,engine);


filterEfficiency(OutEvt.electrons(),ElectronEff,engine);
filterEfficiency(OutEvt.muons(),MuonEff,engine);

smearElectronEnergy(OutEvt.electrons(),engine);
smearMuonMomentum(OutEvt.muons(),engine);
*/


//cout << "After smearing we have " << OutEvt.electrons().size() << " electrons and " << OutEvt.muons().size() << " muons" << endl;
//return OutEvt;




}
