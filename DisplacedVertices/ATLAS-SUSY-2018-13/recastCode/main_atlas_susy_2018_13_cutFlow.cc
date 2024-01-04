// An example recasting code for applying the selection from ATLAS-SUSY-2016-08 and computing the efficiency

#include <iostream>
#include "Pythia8/Pythia.h"
#include "helperFunctions.h"
#include "Pythia8Plugins/HepMC2.h"
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>




using namespace Pythia8;


int run(int nevents, const string & cfgfile, const string & infile, const string &inifile)
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
  Event& event = pythia.event;
  pythia.readFile( cfgfile );
  if ( infile.find(".slha") != std::string::npos ){
    cout << "Using SLHA file as input" << endl;
    pythia.readString("SLHA:file = " + infile);
    if ( nevents < 0) {
    	nevents = 100;
    	cout << "Negative number of events for SLHA input. Setting nevents to " << nevents << endl;
    }
  }
  else{
    cout << "Using LHE file as input" << endl;
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + infile);
  }


  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io("pythia.hepmc", std::ios::out);


  //Get access to particle data:
  ParticleData& pData = pythia.particleData;

  // Initialize.
  pythia.init();

  //Define jet clustering
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rjet);

  int iAbort = 0;
  double nCuts = 0.;
  double nCutsMET = 0.;
  double nCutsMETeff = 0.;
  double nCutsDV = 0.;

  // Begin event loop.
  for (int iEvent = 0; iEvent < nevents; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < 10) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }


    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;


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
    double MET = getMissingMomentum(event).pT();
    if (MET < METcut) {continue;}

    nCutsMET += 1.;

    //Apply pres-selection DV cut:
    if (DVs.size() < 1){continue;}

    nCutsDV += 1.;

    //Get event selection efficiency:
    double Rmax = 0.;
    for (int i =0; i < DVs.size(); ++i){
        Rmax = max(Rmax,DVs[i].vDec().pT());
    }
    double evEff = getEvEff(MET,Rmax);

    nCutsMETeff += evEff;

    //Get probability for reconstructing at least one DV:
    double dvEff = 1.;
    for(int i = 0; i < DVs.size(); ++i){
        dvEff *= 1.-DVs[i].DVeff;
    }
    dvEff = 1. - dvEff;

    //Compute event weight:
    double eventWeight = evEff*dvEff;

    nCuts += eventWeight;


  // End of event loop.
  }

  float evtNorm = 32./float(nevents);
  cout << " Initial number of events = " << float(nevents)*evtNorm << endl;
  cout << " Number of events after pre-selection MET cut = " << nCutsMET*evtNorm << endl;
  cout << " Number of events after pre-selection DV cut = " << nCutsDV*evtNorm << endl;
  cout << " Number of events after MET selection = " << nCutsMETeff*evtNorm << endl;
  cout << " Number of events after DV selection = " << nCuts*evtNorm << endl;
  cout << " Signal Efficiency = " << nCuts/float(nevents)
                           << " ( " << nCuts << " evts )" << endl;

  // Done.
  return 0;
}


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <slhafile/lhefile>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <slhafile/lhefile>:  input SLHA or LHE file [test.slha]" << endl;
	  cout << "        -c <pythia config file>:  pythia config file [pythia8.cfg]" << endl;
	  cout << "        -p <parameters file>:  parameters file [parameters.ini]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  int nevents = 100;
  string slhafile = "test.slha";
  string outfile = "test.eff";
  string cfgfile = "pythia8.cfg";
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
      slhafile = argv[i+1];
      i++;
      continue;
    }



    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(nevents, cfgfile, slhafile, parfile);

  return 0;
}
