// An example recasting code for applying the selection from ATLAS-SUSY-2016-08 and computing the efficiency

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
	double MET = confFile.get<double>("Cuts.MET");
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
  int nCuts = 0;

  // Begin event loop.
  for (int iEvent = 0; iEvent < nevents; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < 10) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    bool passCuts;
    //Split analysis is two bunchs: 75% and 25%
    float lumCut = (std::rand()/(float)RAND_MAX); //lumCut = random(0,1)
    if (lumCut < jetCutsFraction){
        passCuts = applyJetCuts(event, jetDef, pTjet,
            maxJetChargedPT,minJetPt1, minJetPt2,minPVdistance);
    }
    else {passCuts = true;}
    if (!passCuts){continue;}

    //Get good displaced vertices:
    vector<DisplacedVertex> DVs = getDVs(event,minPVdistance,maxRDV,maxZDV,
                                        minTrackPT, minTrackD0, minDecProd,	minDVmass);

    passCuts = applyCuts(event, MET, DVs);

    if (!passCuts) continue;
    nCuts += 1;


  // End of event loop.
  }

  // Final statistics, flavor composition and histogram output.
//  pythia.stat(); config file [pythia8.cfg]
  cout << " Efficiency = " << float(nCuts)/float(nevents)
                           << " ( " << nCuts << " evts )" << endl;
  fprintf(OutputFile,"Efficiency = %1.3e, Total Number of Events = %i, Number of Events after cuts = %i \n",float(nCuts)/float(nevents),nevents,nCuts);
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
