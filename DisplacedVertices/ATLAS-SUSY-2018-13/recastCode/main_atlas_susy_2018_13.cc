// An example recasting code for applying the selection from ATLAS-SUSY-2018-13 and computing the efficiency

#include <iostream>
#include "Pythia8/Pythia.h"
#include "helperFunctions.h"
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>



using namespace Pythia8;


int run(int nevents, const string & cfgfile, const string & slhafile, const string &inifile, const string & outfile)
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
  FILE* OutputFile = fopen(outfile.c_str(), "w");
  // Generator. Shorthand for the event.
  Pythia pythia("",false); //Set printBanner to false
  Event& event = pythia.event;
  pythia.readFile( cfgfile );
  pythia.readString("SLHA:file = " + slhafile);

  //Get access to particle data:
  ParticleData& pData = pythia.particleData;

  // Initialize.
  pythia.init();

  //Define jet clustering
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rjet);

  int iAbort = 0;
  double nCuts = 0.;

  // Begin event loop.
  for (int iEvent = 0; iEvent < nevents; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < 10) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    bool passCuts;

    std::vector<DisplacedVertex> dvs = getDVs(event);
    // All jets
    vector <fastjet::PseudoJet> jets = getJets(event,25.0,5.0);
    // Displaced jet candidates
    vector <fastjet::PseudoJet> jetCandidates = getJets(event,25.0,2.5);
    // Check if any of the jet constituents is a LLP daughter
    


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

    //Apply pres-selection DV cut:
    if (DVs.size() < 1){continue;}

    //Get event selection efficiency:
    double Rmax = 0.;
    for (int i =0; i < DVs.size(); ++i){
        Rmax = max(Rmax,DVs[i].vDec().pT());
    }
    double evEff = getEvEff(MET,Rmax);

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

  cout << " Efficiency = " << nCuts/float(nevents)
                           << " ( " << nCuts << " evts )" << endl;
  fprintf(OutputFile,"Efficiency = %1.3e, Total Number of Events = %i, Effective Number of Events after cuts = %1.3e \n",nCuts/float(nevents),nevents,nCuts);
  fclose(OutputFile);


  // Done.
  return 0;
}


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <slhafile>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <slhafile>:  input SLHA file [test.slha]" << endl;
	  cout << "        -o <output file>:  output file name [test.eff]" << endl;
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

    if ( s== "-o" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      outfile = argv[i+1];
      i++;
      continue;
    }


    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(nevents, cfgfile, slhafile, parfile, outfile);

  return 0;
}
