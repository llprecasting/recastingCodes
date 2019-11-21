// Reads an input LHE or an input SLHA file and generate events using Pythia 8.
// The efficiencies for 1Cand-FullDetector and 2Cand-FullDetector searches are displayed and
// stored in the output file.

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
  // SingleMuTurnOn - provided turn-on histogram for single-muon trigger
  // LooseEff - provided efficiency histogramm for loose candidates
  // TightPromotionEff - provided efficiency histogramm for promoting loose candidates to tight ones
  // MToFFullDet - provided mass resolution histogram
  // decayInsideAtlas - function to check for a decay vertex inside the ATLAS detector
  // lowerLimit_MToF - lower mass limit for signal region (ToF)
  // Particles - list of generator-level particles in the event
  // trandom - e.g. a TRandom3 object


  //Load the ATLAS histograms:
  TFile* InputFile = new TFile("recastCode/ATLAS_data/HEPData-ins1718558-v2-root.root");
  TH1F* EtmissTurnOn = (TH1F*) InputFile->GetDirectory("Table 22")->Get("Hist1D_y1");
  TH2F* SingleMuTurnOn = (TH2F*) InputFile->GetDirectory("Table 23")->Get("Hist2D_y1");
  TH2F* LooseEff = (TH2F*) InputFile->GetDirectory("Table 25")->Get("Hist2D_y1");
  TH2F* TightPromotionEff = (TH2F*) InputFile->GetDirectory("Table 26")->Get("Hist2D_y1");
  TH1F* MToFFullDet = (TH1F*) InputFile->GetDirectory("Table 29")->Get("Hist1D_y1");
  TH1F* MToFFullDetErr = (TH1F*) InputFile->GetDirectory("Table 29")->Get("Hist1D_y1_e1");


  //Create vectors for storing the number of event in each SR/mass window
  std::vector<float> massToFLoose_min  = {175.,375.,600.,825.};
  std::vector<float> massToFTight_min  = {150.,350.,575.,80.};
  std::vector<int> nEvts_SRLoose = {0,0,0,0};
  std::vector<int> nEvts_SRTight = {0,0,0,0};


  // Initialize.
  pythia.init();

  std::mt19937 engine(42);

  int nPass = 0;
  int iAbort = 0;
  // Begin event loop.
  int iEvent = 0;

  int nMET = 0;
  int nCandidate = 0;
  int nMassWindow = 0;
  float avgCandEff = 0.;

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


	std::vector<Particle> candidates;
	int nHSCPs = 0;
	bool TriggerAccept = false;
	// Etmiss trigger accept
	if (Etmiss > 300.) TriggerAccept = true;
	int   bin     = EtmissTurnOn->GetXaxis()->FindBin(Etmiss);
	float eff_Met = EtmissTurnOn->GetBinContent(bin);
	float metRandom = (std::rand()/(float)RAND_MAX);
	if (metRandom < eff_Met) TriggerAccept = true;

    for (int i = 0; i < event.size(); ++i) {

    	Particle particle = event[i];
		// isHSCP identifies the HSCP candidates in the event record
		if (!isHSCP(particle,event)) continue;
		++nHSCPs;

		//		cout << "index = " << i << endl;
		//		cout << "HSCP = " << particle.name() << endl;
		//		cout << "HSCP charge = " << particle.charge() << endl;
		//		cout << "Decay inside ATLAS = " << decayInsideAtlas(particle) << endl;
		//		cout << "Decay vertex = " << particle.vDec() << endl;

		if (decayInsideAtlas(particle)) continue;
		//Only single charged particles are considered
		if (abs(particle.charge()) != 1) continue;

		float eta  = fabs(particle.eta());
		float beta  = particle.pAbs()/particle.e();
		int bin_eta   = SingleMuTurnOn->GetXaxis()->FindBin(eta);
		int bin_beta  = SingleMuTurnOn->GetYaxis()->FindBin(beta);
		float effTrig   = SingleMuTurnOn->GetBinContent(bin_eta, bin_beta);
		float triggerRandom = (std::rand()/(float)RAND_MAX);
		if (triggerRandom < effTrig) TriggerAccept = true;

		// Consider only particles with minimum transverse momentum and beta and maximum eta
		float eta  = fabs(particle.eta());
		float beta  = particle.pAbs()/particle.e();
		if (beta < 0.2) continue;
		if (particle.pT() < 26.) continue;
		if (eta > 2.5) continue;
		candidates.push_back(particle);
    }
    // End of particle loop.

	if (!TriggerAccept) continue;
	if (candidates.size() == 0) continue;
	++nCandidate;

	std::vector<Particle> looseCandidates;
	std::vector<Particle> tightCandidates;
	for (int i = 0; i < candidates.size(); ++i){

		Particle particle = candidates[i];
		float eta = fabs(particle.eta());
		float beta = particle.pAbs()/particle.e();
		// Estimate efficiencies
		int bin_eta = LooseEff->GetXaxis()->FindBin(eta);
		int bin_beta = LooseEff->GetYaxis()->FindBin(beta);
		float effLoose = LooseEff->GetBinContent(bin_eta, bin_beta);
		float effPromotion = TightPromotionEff->GetBinContent(bin_eta, bin_beta);
		float looseRandom = (std::rand()/(float)RAND_MAX);
		float promoteRandom = (std::rand()/(float)RAND_MAX);

		// Momentum cut loose selection
		if (particle.pAbs() < 100.) continue;
		if (particle.pT() < 70.) continue;

		if (looseRandom < effLoose) {
		    looseCandidates.push_back(particle);
		    // Momentum cut tight selection
		    if (particle.pAbs() < 200.) continue;
		    // Sample tight promotion of candidate
		    if (promoteRandom < effPromotion) {
		    	tightCandidates.push_back(particle);
		    }
		  }
	}

	// Final definition of the different signal regions
	if (looseCandidates.size() == 2){
		float mass1 = fabs(looseCandidates[0].m());
		float mass2 = fabs(looseCandidates[1].m());
	    // Events with exactly two loose candidates
	    // Sample the ToF mass for two full-detector candidates
	    int   bin_massToF1        = MToFFullDet->GetXaxis()->FindBin(mass1);
	    float massToF_mean1       = MToFFullDet->GetBinContent(bin_massToF1);
	    bin_massToF1 			  = MToFFullDetErr->GetXaxis()->FindBin(mass1);
	    float massToF_resolution1 = MToFFullDetErr->GetBinError(bin_massToF1);
	    int   bin_massToF2        = MToFFullDet->GetXaxis()->FindBin(mass2);
	    float massToF_mean2       = MToFFullDet->GetBinContent(bin_massToF2);
	    bin_massToF2		      = MToFFullDetErr->GetXaxis()->FindBin(mass2);
	    float massToF_resolution2 = MToFFullDetErr->GetBinError(bin_massToF2);
	    std::normal_distribution<double> gaussToF1(massToF_mean1, massToF_resolution1);
	    std::normal_distribution<double> gaussToF2(massToF_mean2, massToF_resolution2);
	    float massToF1   = gaussToF1(engine);
	    float massToF2   = gaussToF2(engine);

		for (int j = 0; j < nEvts_SRLoose.size(); ++j){
			// Apply final mass requirements for each SR
			if (std::min(massToF1, massToF2)  < massToFLoose_min[j]) continue;
			++nEvts_SRLoose[j];
		}
	}
	else if (tightCandidates.size() == 1) {
		// Events with exactly one tight candiate and not two loose ones
		float mass = fabs(tightCandidates[0].m());
		// Sample the ToF mass for one full-detector candidate
		int   bin_massToF        = MToFFullDet->GetXaxis()->FindBin(mass);
		float massToF_mean       = MToFFullDet->GetBinContent(bin_massToF);
		float massToF_resolution = MToFFullDet->GetBinError(bin_massToF);
		std::normal_distribution<double> gaussToF(massToF_mean, massToF_resolution);
		float massToF = gaussToF(engine);

		for (int j = 0; j < nEvts_SRTight.size(); ++j){
			// Apply final mass requirements for each SR
			if (massToF  < massToFTight_min[j]) continue;
			++nEvts_SRTight[j];
		}
	}

//    cout << "-------------- Event " << iEvent << "----------------" << endl;
  // End of event loop.
  }

//  pythia.stat();

//  cout << "MET trigger eff = " << float(nMET)/float(iEvent) << endl;
//  cout << "Candidate eff = " << float(nCandidate)/float(iEvent) << "  (average candidate eff = " << avgCandEff/float(recHist->GetEntries()) << ")" << endl;
//  cout << "Mass window eff = " << float(nMassWindow)/float(iEvent) << endl;

  //Save histograms to file:


  for (int i = 0; i < nEvts_SRLoose.size(); ++i){
      fprintf(OutputFile,"(mTOF>%3.0f) total Efficiency (Loose): %1.3e\n", massToFLoose_min[i],
      		float(nEvts_SRLoose[i])/float(iEvent));
      fprintf(stdout,"(mTOF>%3.0f) total Efficiency (Loose): %1.3e\n", massToFLoose_min[i],
      		float(nEvts_SRLoose[i])/float(iEvent));
  }
  for (int i = 0; i < nEvts_SRTight.size(); ++i){
      fprintf(OutputFile,"(mTOF>%3.0f) total Efficiency (Tight): %1.3e\n", massToFTight_min[i],
      		float(nEvts_SRTight[i])/float(iEvent));
      fprintf(stdout,"(mTOF>%3.0f) total Efficiency (Tight): %1.3e\n", massToFTight_min[i],
      		float(nEvts_SRTight[i])/float(iEvent));
  }

  //fprintf(OutputFile,"<\\total>\n");
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
  string cfgfile = "pythia8_stau.cfg";
  string outfile = "test.dat";
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
