#include <ctime>

#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/HepMC2.h"
#include "include/gzstream.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
using namespace Pythia8;



bool registerinHCAL(Pythia8::Particle* particle, const Pythia8::Event &event)
{
  if(!particle->isHadron()) {return false;}
  double abseta=fabs(particle->eta());
  Pythia8::Vec4 prodVertex = particle->vProd();
  // idea is that if the particle is produced inside or before the HCAL, we count it here.
  // Will disregard the parent particles in that case ...
  if(abseta > 1.3)
      {
	
	if (prodVertex.pT() > 2.95) return false;
      }
    else if(abseta < 3.0)
      {
	if (fabs(prodVertex.pz()) > 5.5) return false;
      }
    else if(abseta < 5.2)
      {
	if (fabs(prodVertex.pz()) > 12.0) return false;
      }
    else{
      return false;
    }
  
  
  
  if(particle->isFinal())
    {
    return true;
  }
  else
    {
      // Now have to check decays; it needs to have decayed *after leaving HCAL* or otherwise we count the decay products
      Pythia8::Vec4 vdec = particle->vDec();
      if(abseta > 1.3)
      {
	
	if (vdec.pT() > 2.95) return true;
      }
    else if(abseta < 3.0)
      {
	if (fabs(vdec.pz()) > 5.5) return true;
      }
    else if(abseta < 5.2)
      {
	if (fabs(vdec.pz()) > 12.0) return true;
      }
      
      return false;
      
      
    }
  

  
}

/*
bool ischargedTrack(Pythia8::Particle* particle, const Pythia8::Event &event)
{
  // oh good we're keeping charged tracks too

  if(!particle->isCharged()) {return false;}

  if(abseta > 1.3)
      {
	
	if (prodVertex.pT() > 2.95) return false;
      }
    else if(abseta < 3.0)
      {
	if (fabs(prodVertex.z()) > 5.5) return false;
      }
    else if(abseta < 5.2)
      {
	if (fabs(prodVertex.z()) > 12.0) return false;
      }
    else{
      return false;
    }
  
  Pythia8::Vec4 pos = prt->vProd();
  Pythia8::Vec4 dec = prt->vDec();


}
*/


bool isgoodHadron(Pythia8::Particle* particle, const Pythia8::Event &event)
{
  if(!particle->isHadron()) {return false;}
  //  if(!particle->isCharged()) { return false;}

  if(!particle->isFinal())
  {
    // check it's not part of a shower (happens quickly)
        std::vector<int> daughters = particle->daughterListRecursive();
	      for (int i = 0; i < daughters.size(); i++){
		    int idaughter = daughters[i];
		      const Pythia8::Particle* daughterParticle = &event[idaughter];
		    if (daughterParticle->id() == particle->id()) return false;
	      }

        // now check that it doesn't decay before calorimeter, and that it has a reasonable track length
	
  }

  
  if(!registerinHCAL(particle,event)) {return false;}

  // condition to be in the HCAL 
  
  //Pythia8::Vec4 prodVertex = particle->vProd();
  //if( (prodVertex.pT() > 1100.) || (abs(prodVertex.pz()) > 3100.)) {return false;}

  /*

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

        // now check that it doesn't decay before calorimeter, and that it has a reasonable track length
	
  }

  */

  return true;
}


void vecout(const Pythia8::Vec4 vec, GZ::ogzstream &ss)
{
  ss << vec.px() << " " << vec.py() << " " <<  vec.pz() << " " << vec.e() << " ";
}

void writeline(Pythia8::Particle* prt, GZ::ogzstream &ss)
{
  // only visible stuff kept here
  Pythia8::Vec4 pos = prt->vProd();
  Pythia8::Vec4 dec = prt->vDec();
  int pid = prt->id();
  int charge=prt->chargeType();
  if(prt->isHadron())
    {
      if(prt->isCharged())
	{
	  ss << "CH ";
	}
      else
	{
	  ss << "NH ";
	}
    }
  else
    {
      ss << "P ";
    }

  ss << pid << " " << charge << " ";
  vecout(prt->p(),ss);
  if(pos.pAbs() > 0.0001)
    {
      vecout(pos,ss);
    }
  else
    {
      ss << "0.0 0.0 0.0 0.0 ";
    }
  
  if(!(prt->isFinal()))
    {
      vecout(dec,ss);
    }
  else
    {
      ss << "0.0 0.0 0.0 0.0";
    }
  ss<< std::endl;

}


int main() {
  int nEvent  = 20000;
  
  //ogzstream* pileupstream= new ogzstream("minbias.pileup.gz") ;

  GZ::ogzstream pileupstream;
  pileupstream.open("minbias.dat.gz");
  //pileupstream = new ogzstream("minbias.pileup.gz");

  
  /*
  HepMC::IO_GenEvent* hepmcevts;
  HepMC::Pythia8ToHepMC ToHepMC;
  hmcgzstream = new ogzstream("minbias.hepmc.gz");
  hepmcevts = new HepMC::IO_GenEvent(*hmcgzstream);*/
  


  
  //string pdfSet = "LHAPDF6:PDF4LHC15_nlo_asvar";

  string tune="Tune:pp = 8"; // A2 tune, MSTW 2008 LO, needs LHAPDF by default
  Pythia pythia;
  Event& event = pythia.event;
 

  double eCM =  13000.;

  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("SoftQCD:doubleDiffractive = on");
  pythia.readString("Next:numberShowEvent = 1");
  pythia.readString("Tune:preferLHAPDF = 2");
  
  pythia.readString(tune);
  //pythia.readString("PDF:pSet = " + pdfSet);
  pythia.settings.parm("Beams:eCM", eCM);

 
  
  pythia.init();
  int iEvent=0;
  int nAbort=20;
  int iAbort=0;
  while(iEvent < nEvent)
    {

      // Generate events.  Skip if error.
      if (!pythia.next()) {
	iAbort++;
	if(iAbort < nAbort) continue;

	std::cout << "Too many errors, aborting!" << std::endl;
	break;
	  
      }
      iEvent++;
      //pileupstream << "asdasda"<< std::endl;

      // now let's start the storage business

      pileupstream << "E " << iEvent << std::endl;

      for(int npart = 0; npart< event.size(); npart++)
	{
	  Pythia8::Particle* prt = &event[npart];
	  int pid=prt->id();
	  int apid=abs(pid);
	  // check if is final and produced inside tracker or if it decays outside enough (not yet)
	  //cout << pid << endl;
	  if(prt->isFinal())
	    {
	      
	      if(prt->isVisible())
		{

		  if(prt->isHadron())
		    {
		      if(registerinHCAL(prt,event))
			{
			  writeline(prt,pileupstream);
			}
			continue;
		    }
		      else
			{
			  Pythia8::Vec4 prodVertex = prt->vProd();
			  if( (prodVertex.pT() > 7000.) || (abs(prodVertex.pz()) > 12000.)) continue;
		  
		  
			  writeline(prt,pileupstream);

			}
		  
		}
	      
             
             
         }
	  else
	    {
	      // cout << "is not final ";
	      if(isgoodHadron(prt,event))
		{

		  //  if(prt->isCharged())
		  //	{
	 //HEP::Particle* newpart = new HEP::Particle(prt->px(),prt->py(),prt->pz(),prt->e(),pid);
        //Pythia8::Vec4 pos = prt->vProd();
		  // Pythia8::Vec4 vdec = prt->vDec();
		  // if(vdec.pT() < 33.0) continue; // check to make sure it doesn't decay too soon
	writeline(prt,pileupstream);
	//  }

		}
	
    
	    }

	}

 
      
      //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      //ToHepMC.fill_next_event( pythia, hepmcevt );
      // Write the HepMC event to file. Done with it.
      //ascii_io << hepmcevt;
      //hepmcevts << hepmcevt;
      
      //delete hepmcevt;
      
  }
  double sigma = pythia.info.sigmaGen();
  pythia.stat();
   cout << "xs: " << sigma << endl;
  
  return 0;
}
