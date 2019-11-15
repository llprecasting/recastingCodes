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
#include "TH1.h"
#include "TH2.h"
#include <string>
#include <random>



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

  // EtmissTurnOn - provided turn-on histogram for Etmiss trigger
  // IDCaloEff - provided efficiency histogram for MS-agnostic candidates
  // MToFMSagno - provided mass resolution histogram (dE/dx)
  // MdedxMSagno - provided mass resolution histogram (ToF)
  // isPrimaryRHadron - function identifying the initial R-Hadrons in the event record
  // RhadCharge - function returning the charge of an R-hadron
  // decayBeforeEndHcal - function to check for a decay vertex inside or before tile calorimeter
  // lowerLimit_MdEdx - lower mass limit for signal region (dE/dx)
  // lowerLimit_MToF - lower mass limit for signal region (ToF)



  //Load the ATLAS histograms:
  TFile* InputFile = new TFile("recastCode/ATLAS_data/HEPData-ins1718558-v2-root.root");
  TH1F* EtmissTurnOn = (TH1F*) InputFile->GetDirectory("Table 22")->Get("Hist1D_y1");
  TH2F* IDCaloEff = (TH2F*) InputFile->GetDirectory("Table 24")->Get("Hist2D_y1");
  TH1F* MToFMSagno = (TH1F*) InputFile->GetDirectory("Table 28")->Get("Hist1D_y1");
  TH1F* MdedxMSagno = (TH1F*) InputFile->GetDirectory("Table 27")->Get("Hist1D_y1");

  //Create vectors for storing the number of event in each SR/mass window
  std::vector<float> massToF_min  = {350.,550.,700.,850.};
  std::vector<float> massdEdx_min = {300.,450.,600.,750.};
  std::vector<int> nEvts_SR = {0,0,0,0};


  // Initialize.
  pythia.init();

  std::mt19937 engine(42);

  int nPass = 0;
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

    ++iEvent;


    // All events
    double Etmiss = getMissingMomentum(event).pT();
    //cout << "MET = " << Etmiss << endl;
    // Trigger decision
    if (Etmiss < 300.) {
        int   bin     = EtmissTurnOn->GetXaxis()->FindBin(Etmiss);
        float eff_Met = EtmissTurnOn->GetBinContent(bin);
        float metRandom = (std::rand()/(float)RAND_MAX); //lumCut = random(0,1)
        if (metRandom > eff_Met) continue;
//        cout << "Passed MET cut: " << metRandom << " > " << Etmiss << endl;
    }
//    cout << "Passed MET" << endl;

    // Events that passed the trigger
	std::vector<Particle> candidates;
    for (int i = 0; i < event.size(); ++i) {

    	Particle particle = event[i];

		// isPrimaryRHadron identifies the initial R-Hadrons in the event record
		if (!isPrimaryRHadron(particle)) continue;

//		cout << "index = " << i << endl;
//		cout << "R-hadron = " << particle.name() << endl;
//		cout << "R-hadron charge = " << particle.charge() << endl;
//		cout << "Decay before endCal = " << decayBeforeEndHcal(particle) << endl;
//		cout << "Decay vertex = " << particle.vDec() << endl;

		// RhadCharge returns the charge of an R-hadron,
		// so we consider only R-hadrons charged after hadronisation
		if (abs(particle.charge()) != 1) continue;



		// decayBeforeEndHcal checks for a decay vertex before the end of the tile calorimeter,
		// so we only consider R-hadrons that decay after that (see measures above)
		if (decayBeforeEndHcal(particle)) continue;

//		cout << "Good candidate = " << particle.name() << " pT = " << particle.pT() << " p = " << particle.pAbs() << endl;

		// Consider only particles with minimum (transverse) momentum
		if (particle.pAbs() < 200.) continue;
		if (particle.pT() < 50.) continue;

		//Estimated decision
		float eta  = particle.eta();
		float beta     = particle.pAbs()/particle.e();
		int   bin_eta  = IDCaloEff->GetXaxis()->FindBin(fabs(eta));
		int   bin_beta = IDCaloEff->GetYaxis()->FindBin(beta);
		float effCand  = IDCaloEff->GetBinContent(bin_eta, bin_beta);

		float candRandom = (std::rand()/(float)RAND_MAX); //lumCut = random(0,1)
//		cout << "beta = " << beta << " eta = " << eta << endl;
//		cout << "Calorimeter eff = " << effCand << " random = " << candRandom << endl;
		if (candRandom > effCand) continue;

		candidates.push_back(particle);
	// End of particle loop.
    }

	if (candidates.size() == 0) continue;

	float Mass = 0.0;
	std::vector<int> passed_SR = {0,0,0,0};
	for (int i = 0; i < candidates.size(); ++i){

		Mass = candidates[i].m();
		// Events with at least one candidate
		// Sample the masses and apply the final mass window requirements
		int bin_massToF  = MToFMSagno->GetXaxis()->FindBin(Mass);
		int bin_massdEdx = MdedxMSagno->GetXaxis()->FindBin(Mass);

		// Sample the ToF mass for ID+Calo candidates
		float massToF_mean = MToFMSagno->GetBinContent(bin_massToF);
		float massToF_resolution = MToFMSagno->GetBinError(bin_massToF);
		std::normal_distribution<double> gaussToF(massToF_mean, massToF_resolution);
		float massToF = gaussToF(engine);

		 // Sample the dEdx mass for ID+Calo candidates
		float massdEdx_mean       = MdedxMSagno->GetBinContent(bin_massdEdx);
		float massdEdx_resolution = MdedxMSagno->GetBinError(bin_massdEdx);
		std::normal_distribution<double> gaussdEdx(massdEdx_mean, massdEdx_resolution);
		float massdEdx   = gaussdEdx(engine);

		for (int j = 0; j < nEvts_SR.size(); ++j){
			// Apply final mass requirements for each SR
			if (massToF < massToF_min[j] || massdEdx < massdEdx_min[j]) continue;
			++passed_SR[j];
		}
	}
	// If passed_SR, add event to signal region:
	for (int j = 0; j < nEvts_SR.size(); ++j){
		// Apply final mass requirements for each SR
		if (passed_SR[j] == 0) continue;
		++nEvts_SR[j];
	}


//    cout << "-------------- Event " << iEvent << "----------------" << endl;
  // End of event loop.
  }

  pythia.stat();

  for (int i = 0; i < nEvts_SR.size(); ++i){
        //fprintf(OutputFile,"M>%3i Total Efficiency: %1.3e +- %1.3e\n",Mi*100,totalEffs[Mi].first/iEvent,err);
        fprintf(stdout,"(mTOF>%3.0f, mdEdx > %3.0f) total Efficiency: %1.3e\n", massToF_min[i],massdEdx_min[i],
        		float(nEvts_SR[i])/float(iEvent));
  }
  //fprintf(OutputFile,"<\\total>\n");
  //fclose(OutputFile);
    

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
