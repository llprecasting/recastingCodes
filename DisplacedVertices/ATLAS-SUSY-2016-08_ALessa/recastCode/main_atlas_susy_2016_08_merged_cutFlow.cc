// An example recasting code for applying the selection from ATLAS-SUSY-2016-08 and computing the efficiency
// It should be used with an LHE file containing p p > go go and p p > go go j events
// The MLM matching algorithm is used according to the parameters in pythia8_MLM.cfg
// The MET spectrum (before any cuts) and a python code showing how to plot it are also saved in a .dat file.

#include <iostream>
#include "Pythia8/Pythia.h"
#include "helperFunctions.h"
#include "Pythia8Plugins/HepMC2.h"
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "Pythia8Plugins/CombineMatchingInput.h"



using namespace Pythia8;


int run(int nevents, const string & cfgfile, const string & lhefile, const string &inifile)
{

  //Get cuts and options from config file:
  boost::property_tree::ptree confFile;
  boost::property_tree::ini_parser::read_ini(inifile, confFile);
  //Base selection cuts:
  double minPVdistance = confFile.get<double>("BaseSelection.minPVdistance");
  double maxRDV = confFile.get<double>("BaseSelection.maxRDV");
  double maxZDV = confFile.get<double>("BaseSelection.maxZDV");
  double pTjet = confFile.get<double>("BaseSelection.pTjet");
  double Rjet = confFile.get<double>("BaseSelection.Rjet");
  //SR cuts:
  double minTrackPT = confFile.get<double>("Cuts.minTrackPT");
  double minTrackD0 = confFile.get<double>("Cuts.minTrackD0");
  double minDVmass = confFile.get<double>("Cuts.minDVmass");
  int minDecProd = confFile.get<int>("Cuts.minDecProd");
  double METcut = confFile.get<double>("Cuts.MET");
  double maxJetChargedPT = confFile.get<double>("Cuts.maxJetChargedPT");
  double minJetPt1 =  confFile.get<double>("Cuts.minJetPt1");
  double minJetPt2 =  confFile.get<double>("Cuts.minJetPt2");
  //Options:
  double jetCutsFraction = confFile.get<double>("Options.applyJetCuts");

  std::srand(500);
  // Generator. Shorthand for the event.
  Pythia pythia("",false); //Set printBanner to false
  pythia.readFile(cfgfile);
  Event& event = pythia.event;

  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhefile);

  // Create UserHooks pointer. Stop if it failed. Pass pointer to Pythia.
  CombineMatchingInput combined;
  UserHooks* matching = combined.getHook(pythia);
  if (!matching) return 1;
  pythia.setUserHooksPtr(matching);

  //Get access to particle data:
  ParticleData& pData = pythia.particleData;

  //Define jet clustering
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rjet);

  //Book Histogram
  Hist hMET("MET", 40, 0., 1000.);

  // Cross section an error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;
  double nCount = 0.;
  double nCuts = 0.;
  double nCutsMET = 0.;
  double nCutsMETeff = 0.;
  double nCutsDV = 0.;

  int nAccepted = 0;

  pythia.init();

  for (int iEvent = 0; iEvent < nevents;  ++iEvent) {

    // Generate events. Quit if at end of file or many failures.
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) {
        cout << "Info: end of input file reached" << endl;
        break;
      }
      else continue;
    }

    //Event weight
    double weight = 1.;

    nCount += weight;
    ++nAccepted;

    double MET = getMissingMomentum(event).pT();
    hMET.fill(MET);

    bool passCuts;
    //Split analysis is two bunchs: 75% and 25%
    float lumCut = (std::rand()/(float)RAND_MAX); //lumCut = random(0,1)
    if (lumCut < jetCutsFraction){
        passCuts = applyJetCuts(event, jetDef, pTjet,
            maxJetChargedPT,minJetPt1, minJetPt2,minPVdistance);
    }
    else {passCuts = true;}
    if (!passCuts){continue;}

    //Get displaced vertex candidates:
    vector<DisplacedVertex> DV_candidates = getDVs(event,minTrackPT, minTrackD0);

    //Get good DVs:
    vector<DisplacedVertex> DVs;
    DisplacedVertex DV;
    Vec4 vertex;
    for (int i=0; i < DV_candidates.size(); ++i){
        DV  = DV_candidates[i];
        vertex = DV.vDec();
        //Apply basic selection efficiency
        if (vertex.pT() < minPVdistance) continue;  //Transverse plane separation from PV > 4mm
        if (vertex.pT() > maxRDV) continue;  //|R_DV| < 300 mm
        if (fabs(vertex.pz()) > maxZDV) continue;  //|z_DV| < 300 mm

        //Apply vertex cuts:
        if (DV.decayProducts.size() < minDecProd){continue;}
        if (DV.mDV() < minDVmass){continue;}

        //Get DV reconstruction efficiency:
        DV.DVeff = getDVEff(DV.mDV(),DV.decayProducts.size(),vertex.pT());
        DVs.push_back(DV);
    }

    //Apply pre-selection MET cut:
    if (MET < METcut) {continue;}

    nCutsMET += 1.*weight;

    //Apply pres-selection DV cut:
    if (DVs.size() < 1){continue;}

    nCutsDV += 1.*weight;

    //Get event selection efficiency:
    double Rmax = 0.;
    for (int i =0; i < DVs.size(); ++i){
        Rmax = max(Rmax,DVs[i].vDec().pT());
    }
    double evEff = getEvEff(MET,Rmax);

    nCutsMETeff += evEff*weight;

    //Get probability for reconstructing at least one DV:
    double dvEff = 1.;
    for(int i = 0; i < DVs.size(); ++i){
        dvEff *= 1.-DVs[i].DVeff;
    }
    dvEff = 1. - dvEff;

    //Compute event weight:
    nCuts += evEff*dvEff*weight;
  // End of event loop.
  }

  pythia.stat();
  delete matching;

  hMET *= 1./float(hMET.getEntries()); //Normalize histogram
  HistPlot hpl("METplot");
  hpl.frame( "MET_Histogram");
  hpl.add(hMET, "h", "$E_{T}^{miss}$");
  hpl.plot();

  cout << "Accepted events = " << nAccepted << " Total xsec after matching (pb) = " << pythia.info.sigmaGen()*1e9 << endl;

  float evtNorm = 32./nCount;
  cout << " Initial number of events = " << nCount*evtNorm << endl;
  cout << " Number of events after pre-selection MET cut = " << nCutsMET*evtNorm << endl;
  cout << " Number of events after pre-selection DV cut = " << nCutsDV*evtNorm << endl;
  cout << " Number of events after MET selection = " << nCutsMETeff*evtNorm << endl;
  cout << " Number of events after DV selection = " << nCuts*evtNorm << endl;
  cout << " Signal Efficiency = " << nCuts/nCount
                           << " ( " << nCuts << " evts )" << endl;

  // Done.
  return 0;
}



void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <slhafile>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <lhefileRoot>: (root) name of the LHE file [test]" << endl;
	  cout << "        -c <pythia config file>:  pythia config file [pythia8.cfg]" << endl;
	  cout << "        -p <parameters file>:  parameters file [parameters.ini]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  int nevents = 100;
  int njet = 1;
  string lhefile = "test";
  string outfile = "test.eff";
  string cfgfile = "pythia8_MLM.cfg";
  string parfile = "parameters.ini";
  for ( int i=1; i!=argc ; ++i )
  {
    string s = argv[i];
    if ( s== "-h" )
    {
      help ( argv[0] );
    }

    if ( s== "-c" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      cfgfile = argv[i+1];
      i++;
      continue;
    }

    if ( s== "-p" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      parfile = argv[i+1];
      i++;
      continue;
    }

    if ( s== "-n" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      nevents = atoi(argv[i+1]);
      i++;
      continue;
    }


    if ( s== "-f" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      lhefile = argv[i+1];
      i++;
      continue;
    }



    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(nevents, cfgfile, lhefile, parfile);

  return 0;
}
