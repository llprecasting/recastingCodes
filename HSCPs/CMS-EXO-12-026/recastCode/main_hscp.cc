// Reads an input LHE or an input SLHA file and generate events using Pythia 8.
// The HSCP efficiencies for each event as well as the 4-momentum of each isolated HSCP
// are stored in the output in a simplified LHE format.
// All the HSCPs are assumed to stable.
// For including finite lifetime effects, the events must be reweighted by Flong = exp(-width*l_out/gamma*beta).

#include <iostream>
#include "Pythia8/Pythia.h"
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include "helperFunctionsHSCP.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH3.h"
#include <string>



using namespace Pythia8;



int run(const string & infile, int nevents, const string & cfgfile, const string & outfile, double width = 0.)
{


  std::srand(500);
  FILE* OutputFile = fopen(outfile.c_str(), "w");
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
  //Load the CMS histograms:
  int nM = 4;
  TFile* InputFile = new TFile("recastCode/CMS_data/EXO13006_ASCII_TABLES_ALL.root");
  TH3F* histoTrigger = (TH3F*)InputFile->Get("P_On(k)");
  TH3F** histoOffline = new TH3F*[nM];
  for(int Mi = 0; Mi < nM; Mi++){
      char histoName[256];
      sprintf(histoName, "P_Off(k,%iGeV)", Mi*100);
      histoOffline[Mi] = (TH3F*)InputFile->Get(histoName);
  };

//  pythia.readString("1000024:m0 = 500.");
  // Initialize.
  pythia.init();

  std::vector< std::pair<double,double> > totalEffs(4,std::make_pair(0.,0.));


  int iAbort = 0;
  // Begin event loop.
  int iEvent = 0;
  while (iEvent < nevents or nevents < 0){

    // If failure because reached end of file then exit event loop.
    if (pythia.info.atEndOfFile()) break;

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < 10) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    //Get isolated HSCPs
    std::vector<Particle*> IsoHSCPs = getIsolatedHSCPs(pythia.event);
	fprintf(OutputFile,"<event>\niev: %d\nnHSCPs: %lu\n",iEvent,IsoHSCPs.size());
    std::vector< std::vector< std::pair<double,double> > > allEffs;
    std::vector<double> Flong(IsoHSCPs.size(),0.);
	for (int i = 0; i < IsoHSCPs.size(); ++i) {
        Particle* HSCP = IsoHSCPs[i];
        std::vector< std::pair<double,double> > HSCPeff = computeHSCPeff(HSCP,histoTrigger,histoOffline);
        Flong[i] = computeFlong(HSCP,width);
        allEffs.push_back(HSCPeff);
        //Print HSCP ID and momentum
        fprintf(OutputFile,"%9d %1.3e %1.3e %1.3e %1.3e %1.3e ",HSCP->id(),HSCP->px(),HSCP->py(),HSCP->pz(),HSCP->e(),HSCP->m0());
        //Print HSCP efficiencies
        for (int j = 0; j < HSCPeff.size(); ++j){
            fprintf(OutputFile,"%1.3e +- %1.3e    ",HSCPeff[j].first,HSCPeff[j].second);
        }
        fprintf(OutputFile,"\n");
    }

    std::vector< std::pair<double,double> > effs = computeTotalEff(allEffs,Flong);
    for (int Mi = 0; Mi < effs.size(); ++Mi){
        totalEffs[Mi].first += effs[Mi].first;
        totalEffs[Mi].second += effs[Mi].second;
    }
    fprintf(OutputFile,"<\\event>\n");
    ++iEvent;

  // End of event loop.
  }

  pythia.stat();

  fprintf(OutputFile,"<total>\n");
  for (int Mi = 0; Mi < totalEffs.size(); ++Mi){
        double err = sqrt(pow(sqrt(totalEffs[Mi].second)/iEvent,2) + pow(totalEffs[Mi].first*sqrt(iEvent)/pow(iEvent,2),2));
        fprintf(OutputFile,"M>%3i Total Efficiency: %1.3e +- %1.3e\n",Mi*100,totalEffs[Mi].first/iEvent,err);
        fprintf(stdout,"M>%3i Total Efficiency: %1.3e +- %1.3e\n",Mi*100,totalEffs[Mi].first/iEvent,err);
  }
  fprintf(OutputFile,"<\\total>\n");
  fclose(OutputFile);
    

  // Done.
  return 0;
}


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <input file>] [-o <output file>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <input file>:  pythia input LHE or SLHA file" << endl;
	  cout << "        -c <pythia config file>:  pythia config file [pythia8.cfg]" << endl;
	  cout << "        -o <output file>:  pythia output LHE file [test.lhe]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]. If n < 0, it will run over all events in the LHE file" << endl;
	  cout << "        -w <width (Gev)>: Optional width to be used [0]" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  float weight = 1.;
  int nevents = 100;
  double width = 0.;
  string cfgfile = "pythia8.cfg";
  string outfile = "test.lhe";
  string infile = "";
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


    if ( s== "-n" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      nevents = atoi(argv[i+1]);
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


    if ( s== "-f" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      infile = argv[i+1];
      i++;
      continue;
    }


    if ( s== "-w" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      width = atof(argv[i+1]);
      i++;
      continue;
    }


    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(infile, nevents, cfgfile, outfile, width);

  return 0;
}
