//////////////////////////////////////////////////////
// ATLAS Displaced Vertex 13 TeV recast WITH TRUTH DV
// Written by Giovanna Cottin (gfcottin@gmail.com)
///////////////////////////////////////////////////
// ATLAS DV+MET displaced vertex analysis in SUSY-2016-08
// Cuts: 
// Event Selection :
//   *   MET > 250 GeV
// DV Selection - Event pass if at least one displaced vertex with :
//   *   pT of tracks coming from DV  > 1 GeV
//   *   |d0| of tracks coming from DV  > 2 mm
//   This last two cuts define the signal region :
//   *   Number of tracks coming out of vertex > 5 
//   *   Invariant mass of vertex > 10 GeV 

#include "Pythia8/Pythia.h"
// Toy Detector Simulation provided by Nishita Desai, needed for jet cuts and MET reco
// #include "ToyDetector-ATLAS.h"
// Toy Detector Simulation modyfied by Giovanna Cottin for tracklessjets
#include "ToyDetector-ATLAS-tracklessjet.h"// both jets and trackless jets
// NEW ATLAS parametrized_truthEff 
#include "parametrized_truthEff.h"

#include <vector>
#include <string>
#include <iterator>
#include <algorithm>

using namespace Pythia8;

int main() {
  Pythia pythia;
  Event& event = pythia.event;
  // Split SUSY model
  pythia.readFile("pythia8.cfg");

 // string ctauValues[12] = {"8.99376000e-01",   "8.99376000e-01", "2.99792000e+00" ,"1.19916800e+01",
 //         "2.99792000e+01",   "1.19916800e+02",   "2.99792000e+02",
 //         "8.99376000e+02",   "2.99792000e+03",   "7.35629610e+03",
 //         "8.99376000e+03",   "1.49896000e+04"};

// For neutralino 1, mg=1400
string massValues[13]={ "100.","100.","200.","400.","600.","800.","1000.","1200.","1270.","1300.","1320.", "1350.", "1370."};

// For neutralino 2, mg=2000
//string massValues[13]={"1.000000e+02","1.000000e+02","2.000000e+02","4.000000e+02","6.000000e+02","8.000000e+02","1.000000e+03","1.200000e+03","1.400000e+03","1.600000e+03","1.800000e+03","1.870000e+03","1.900000e+03"};



  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  ////// Initialize.
  ////// pythia.init();
  //Files
  ofstream effCutFlow;
  effCutFlow.open("Plots/truthDV_data/effCutFlow_1400_varyMass.dat"); 
  ofstream effVsTau;
  effVsTau.open("Plots/truthDV_data/effVsTau_1400_varyMass.dat"); 
  cout<<"###################"<<endl;
  pythia.readString("1000021:m0=1400"); // gluino mass
  cout<<"pythia.particleData.m0(1000021)   = "<<pythia.particleData.m0(1000021)<<endl;
  // Now lifetime is fixed to 1 ns
  pythia.readString("1000021:tau0=299.792");
  cout<<"pythia.particleData.tau0(1000021) = "<<pythia.particleData.tau0(1000021)<<endl;

  /////////////////////
  //Begin mass loop
  int iTau = 0;
  bool isInit = false;
  for (iTau=12; iTau > 0; iTau--){

     // Initialize.
     if(!isInit) {
       pythia.init();
       isInit = true;
     }
     // Consistenty lifetime check
    // width = 3.*1e11*6.581e-25/ctau
     cout<<"###################"<<endl;
     // The gluino is the LLP in this model
    // pythia.readString("1000021:tau0 =" + ctauValues[iTau]);
    pythia.readString("1000022:m0 =" + massValues[iTau]);

    cout<<"pythia.particleData.m0(1000022) = "<<pythia.particleData.m0(1000022)<<endl;

    int iAbort = 0;
    int nPromptMET=0;
    int nJetCuts=0;
    int nDVReco=0;
    int nMaterial=0;
    int nFidutial=0;
    int nTrk=0;
    int nMass=0;
    int nEvtEff=0;
    int nDVEff=0;

     // cout<<"########## BEFORE CALLING DETECTOR #########"<<endl;
    ToyDetector detector;
     // cout<<"########## AFTER CALLING DETECTOR #########"<<endl;

      //Start event loop
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        // cout<<" ####### Event "<< iEvent <<" #########"<<endl;
        // Generate events. Quit if many failure.
        if (!pythia.next()) {
          event.list();
          cout<<" inside pythia!next"<<endl;
          if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
          break;
        }
        if(!detector.getObjects(event)) {
          cout << "No objects found" << endl;
        } 
        else {
           // detector.printObjects();
        }
      bool passEvtPromptMET = false;
      bool passEvtJetCuts   = false; 
      bool passEvtDVReco    = false;
      /////////////////
      //Event Selection
      /////////////////
      // Missing transverse momentum
      // from truth 
      //////////////////////////////
      Vec4 truthMET;
      for (int i = 0; i < event.size(); ++i) {
        if (!event[i].isFinal()) continue;
        if(event[i].idAbs() == 12 || event[i].idAbs() == 14 || event[i].idAbs() == 16
          || event[i].idAbs() == 1000022 || event[i].idAbs() == 1000023 || event[i].idAbs() == 1000025
          || event[i].idAbs() == 1000035){
          truthMET += event[i].p();
        }
      }
      if(truthMET.pT()>200.){
       passEvtPromptMET=true;
      }
      if(passEvtPromptMET) nPromptMET++;
       else continue;

      std::vector<double> JetpT;
      // Number of jets in event
      // cout<<"detector.jets.size()          = "<<detector.jets.size()<<endl;
      // cout<<"detector.tracklessjets.size() = "<<detector.tracklessjets.size()<<endl;

      // Trackless jet cuts
      // Cuts for jets to are applyied only 75 % of the events
      if(iEvent<=0.75*nEvent){
        for (unsigned long j=0; j < detector.tracklessjets.size(); j++) {
          JetpT.push_back(detector.tracklessjets.at(j).pT());
          if(abs(detector.tracklessjets.at(j).eta())<2.8) continue; 
        }
        //Here I cut on the corresponding n jets and pT (ordered in pT)
        bool is1j70AccepCut = (JetpT.size() > 0) && (JetpT[0] > 70.);
        bool is2j25AccepCut = (JetpT.size() > 1) && (JetpT[1] > 25.);
        if(is1j70AccepCut || is2j25AccepCut) passEvtJetCuts=true;
         passEvtJetCuts=true;
        if(passEvtJetCuts) nJetCuts++;
          else continue;
      }
      else{
        passEvtJetCuts=true;
      }
      /////////////////////////////////////////////
      //Displaced Vertex truth identification
      //Find all particles coming from a Rhadron
      /////////////////////////////////////////////  
      std::vector<int> motherIndices;   
      //Get Rhadrons from event
      for (int i= 0; i < event.size(); i++){
        if (abs(event[i].status()) == 104){
          motherIndices.push_back(i);
        }
      }
      //At least 1 "reco" DV in event (we expect two per event at truth level)
      if(motherIndices.size()>=1) passEvtDVReco = true;
      if(passEvtDVReco) nDVReco++;
        else continue;
      // Cuts
      bool vtxPassFidutial = false;
      bool vtxPassesNtrk   = false;
      bool vtxPassesMass   = false;
      bool passDVEff       = false;
      bool passEvtEff      = false; 
      bool passDVAll       = false;
      double Rdecay_largest = 0.0;
      // Loop over truth DVs
      for (int i = 0; i < motherIndices.size(); ++i) {
        int DVindex = motherIndices[i];
        double xDV=event[DVindex].xDec() - event[DVindex].xProd();  
        double yDV=event[DVindex].yDec() - event[DVindex].yProd();  
        double zDV=event[DVindex].zDec() - event[DVindex].zProd();  
        double rDV=sqrt(xDV*xDV+yDV*yDV);
        Vec4 total4p;
        Vec4 total4pPion;  
        //Get all daughters
        vector<int> daughters= event[DVindex].daughterListRecursive();// Need this method to get all daughters
        std::vector<int> daughterIndices;
        for (int j = 0; j < daughters.size(); ++j){
          if (!event[daughters[j]].isFinal()) continue;
          if (std::find(daughterIndices.begin(),daughterIndices.end(),daughters[j]) != daughterIndices.end()) continue;
          daughterIndices.push_back(daughters[j]);
        }
        std::vector<int> chargedDaughters;
        //track loop
        for(int trk=0;trk<daughterIndices.size();trk++){
          int trackIndex = daughterIndices[trk];
          bool passTrackQuality  = true; // true unless probe otherwise
          double rTrk=sqrt(pow2((event[trackIndex].vProd()).px())+ pow2((event[trackIndex].vProd()).py()));
          double zTrk=abs((event[trackIndex].vProd()).pz());   
          double phixy = event[trackIndex].vProd().phi();
          double deltaphi = phixy - event[trackIndex].phi();
          if (abs(deltaphi) > 3.1415) deltaphi = 2 * 3.1415 - abs(deltaphi);
          //d0 in the abscence of magnetic field
          //Make sure tracks and DV is displaced
          double d0 = rTrk*sin(deltaphi);
          //select good daughter tracks 
          if(!event[trackIndex].isCharged()) continue;
          if(abs(d0)<2.0 || event[trackIndex].pT()<1.0) continue;
          chargedDaughters.push_back(trackIndex);
          // DV invariant mass calculation with pionMass hypothesis  
          double particlepx = event[trackIndex].px();
          double particlepy = event[trackIndex].py();
          double particlepz = event[trackIndex].pz();
          double particleE = event[trackIndex].e();
          double particleEPion = sqrt(pow2(0.1395700)+pow2(particlepx)+pow2(particlepy)+pow2(particlepz));
          total4p+=Vec4(particlepx,particlepy,particlepz,particleE);
          total4pPion+=Vec4(particlepx,particlepy,particlepz,particleEPion);
        }//close track loop
        double totalmPion=total4pPion.mCalc();
        // Apply DV cuts. AT LEAST ONE DV satisfying all criteria
        vtxPassFidutial = (4.<rDV) && (rDV<300) && (abs(zDV)<300);
        vtxPassesNtrk   = chargedDaughters.size()>=5;
        vtxPassesMass   = totalmPion>10.;
        // Apply DV efficiency
        float iRndNumber1 = (std::rand()/(float)RAND_MAX); 
        passDVEff         = vertexEff_Regions(rDV, totalmPion, chargedDaughters.size())>iRndNumber1;
        if(vtxPassFidutial && vtxPassesNtrk && vtxPassesMass && passDVEff) {
          //Find max of the one/two DVs passing all selection
          if(rDV>Rdecay_largest) {
            Rdecay_largest=rDV;
          }
          passDVAll=true;
          break;// leave vertex for for, as I have found at least one satisfying all criteria        
        }
      }//close vertex loop
      if(passDVAll){
        // cout<<"AT LEAST ONE VERTEX PASSES ALL CONDITIONS"<<endl;
      }
      //Apply event level efficiency using the maximum rDV present in the event.
      float iRndNumber2 = (std::rand()/(float)RAND_MAX); 
      if(eventEff_MET(Rdecay_largest,truthMET.pT())>iRndNumber2) passEvtEff=true;
      //CutFlow count
      if(vtxPassFidutial) nFidutial++;
        else continue;
      if(vtxPassesNtrk) nTrk++;
        else continue;
      if(vtxPassesMass) nMass++;
        else continue;        
      if(passDVEff) nDVEff++;
        else continue;    
      if(passEvtEff) nEvtEff++;
        else continue;        
    }//End of event loop.
    ////////////////////////
    // CutFlow default ATLAS
    ////////////////////////
    // Number of events - Relative Efficiency - Overall efficiency
    effCutFlow <<" CutFlow for ctau, tau = "<< pythia.particleData.tau0(1000021)<<" , "<<pythia.particleData.tau0(1000021)/(299.792) <<endl;
    effCutFlow <<" CutFlow for tau = " << pythia.particleData.tau0(1000021)<<endl;
    effCutFlow <<" All events "        << nEvent     << " " << nEvent*100./nEvent       << " " << nEvent*100./nEvent    << endl;
    effCutFlow <<" MET Cut  "          << nPromptMET << " " << nPromptMET*100./nEvent   << " " << nPromptMET*100./nEvent<< endl;
    effCutFlow <<" Jet Cut  "          << nJetCuts   << " " << nJetCuts*100./nPromptMET << " " << nJetCuts*100./nEvent  << endl;
    effCutFlow <<" DV Reconstruction " << nDVReco    << " " << nDVReco*100./nJetCuts    << " " << nDVReco*100./nEvent   << endl;
    effCutFlow <<" DV Fidutial "       << nFidutial  << " " << nFidutial*100./nDVReco   << " " << nFidutial*100./nEvent << endl;
    effCutFlow <<" DV Ntrk "           << nTrk       << " " << nTrk*100./nFidutial      << " " << nTrk*100./nEvent      << endl;
    effCutFlow <<" DV Mass "           << nMass      << " " << nMass*100./nTrk          << " " << nMass*100./nEvent     << endl;
    effCutFlow <<" DV Eff "            << nDVEff     << " " << nDVEff*100./nMass        << " " << nDVEff*100./nEvent    << endl;
    effCutFlow <<" Evt Eff "           << nEvtEff    << " " << nDVEff*100./nEvtEff      << " " << nEvtEff*100./nEvent   << endl;
    ////////////////////
    //Efficiency Vs Ctau
    ////////////////////
    // effVsTau<<pythia.particleData.tau0(1000021)<<" "<<setprecision(5)<<fixed<<nEvtEff<<" "<<setprecision(5)<<nEvtEff*1.0/nEvent<<endl;
    effVsTau<<pythia.particleData.m0(1000022)<<" "<<setprecision(5)<<fixed<<nEvtEff<<" "<<setprecision(5)<<nEvtEff*1.0/nEvent<<endl;
        // pythia.stat();
  }//close ctau for
  effVsTau.close();
  effCutFlow.close();  

  return 0;
}

