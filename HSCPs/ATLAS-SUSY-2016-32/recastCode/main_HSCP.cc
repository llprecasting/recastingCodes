// Reads an input LHE or an input SLHA file and generate events using Pythia 8.
// The efficiencies for 1Cand-FullDetector and 2Cand-FullDetector searches are displayed and
// stored in the output file.

#include <iostream>
#include "Pythia8/Pythia.h"
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <string>
#include <random>
#include "helperFunctions.h"



using namespace Pythia8;



int run(const string & infile, int nevents, const string & cfgfile, const string & outputfile,  const string & histofile)
{

  //Set output file names
  std::srand(500);
  string outname, histname;
  if (outputfile == ""){
	  size_t lastindex = infile.find_last_of(".");
	  outname = infile.substr(0, lastindex)+".out";
  }
  else{outname = outputfile;}
  if (histofile == ""){
	  size_t lastindex = infile.find_last_of(".");
	  histname = infile.substr(0, lastindex)+".dat";
  }
  else{histname = histofile;}

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


  // Book histograms.
  TH2F *recHist = new TH2F("Reconstruction","", 20, 0., 2.,20,0.,1.);
  TH1F *massLooseHist = new TH1F("mToFLoose","", 100, 0., 1000.);
  TH1F *massTightHist = new TH1F("mToFTight","", 100, 0., 1000.);
  TH1F *pTHist = new TH1F("pT","", 100, 0., 2000.);
  TH1F *betaHist = new TH1F("beta","", 100, 0., 1.);
  TH1F *metHist = new TH1F("MET","", 100, 0., 1000.);


  //Create vectors for storing the number of event in each SR/mass window
  std::vector<float> massToF2Cand_min  = {150.,350.,575.,800.};
  std::vector<float> massToF1Cand_min  = {175.,375.,600.,825.};
  std::vector<int> nEvts_SR1Cand = {0,0,0,0};
  std::vector<int> nEvts_SR2Cand = {0,0,0,0};


  // Initialize.
  pythia.init();

  std::mt19937 engine(42);

  int nPass = 0;
  int iAbort = 0;
  // Begin event loop.
  int iEvent = 0;

  int nTrigger = 0;
  int nCandidate = 0;
  int nMassWindow = 0;
  int nDecay = 0;
  float avgCandEff = 0.;
  int nHSCPs = 0;

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

    metHist->Fill(Etmiss);

	bool TriggerAccept = false;
	// Etmiss trigger accept
	if (Etmiss > 300.) TriggerAccept = true;
	int   bin     = EtmissTurnOn->GetXaxis()->FindBin(Etmiss);
	float eff_Met = EtmissTurnOn->GetBinContent(bin);
	float metRandom = (std::rand()/(float)RAND_MAX);
	if (metRandom < eff_Met) TriggerAccept = true;

	//If the MET trigger was not activated, check for the muon trigger:
	if (!TriggerAccept){
		for (int i = 0; i < event.size(); ++i) {

			Particle particle = event[i];
			// isHSCP identifies the HSCP candidates in the event record
			if (!isHSCP(particle,event)) continue;

			if (decayInsideAtlas(particle)) continue;
			//Only single charged particles are considered
			if (abs(particle.charge()) != 1) continue;

			float eta  = fabs(particle.eta());
			float beta  = particle.pAbs()/particle.e();
			int bin_eta   = SingleMuTurnOn->GetXaxis()->FindBin(eta);
			int bin_beta  = SingleMuTurnOn->GetYaxis()->FindBin(beta);
			float effTrig   = SingleMuTurnOn->GetBinContent(bin_eta, bin_beta);
			float triggerRandom = (std::rand()/(float)RAND_MAX);
			if (triggerRandom < effTrig) {
				TriggerAccept = true;
				break;
			}
		}
	}
	if (!TriggerAccept) continue;
	++nTrigger;

	//Get HSCP candidates:
	std::vector<Particle> looseCandidates;
	std::vector<Particle> tightCandidates;
	for (int i = 0; i < event.size(); ++i) {

		Particle particle = event[i];
		// isHSCP identifies the HSCP candidates in the event record
		if (!isHSCP(particle,event)) continue;
		//Only single charged particles are considered
		if (abs(particle.charge()) != 1) continue;
		++nHSCPs;
		if (decayInsideAtlas(particle)) continue;
		++nDecay;

		float eta  = fabs(particle.eta());
		float beta  = particle.pAbs()/particle.e();
		float pT = particle.pT();
		float p = particle.pAbs();

		//Pre-selection:
		if (eta > 2.0) continue;
		if (pT < 70.0) continue;
		if (p < 100.0) continue;

		betaHist->Fill(beta);
		recHist->Fill(eta,beta);
		pTHist->Fill(pT);

		// Estimate efficiencies
		int bin_eta = LooseEff->GetXaxis()->FindBin(eta);
		int bin_beta = LooseEff->GetYaxis()->FindBin(beta);
		float effLoose = LooseEff->GetBinContent(bin_eta, bin_beta);
		float effPromotion = TightPromotionEff->GetBinContent(bin_eta, bin_beta);
		float looseRandom = (std::rand()/(float)RAND_MAX);
		float promoteRandom = (std::rand()/(float)RAND_MAX);

		avgCandEff += effLoose;

		if (looseRandom < effLoose) {
		    looseCandidates.push_back(particle);
		    ++nCandidate;
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
	    float massToF_resolution1 = MToFFullDetErr->GetBinContent(bin_massToF1);
	    int   bin_massToF2        = MToFFullDet->GetXaxis()->FindBin(mass2);
	    float massToF_mean2       = MToFFullDet->GetBinContent(bin_massToF2);
	    bin_massToF2		      = MToFFullDetErr->GetXaxis()->FindBin(mass2);
	    float massToF_resolution2 = MToFFullDetErr->GetBinContent(bin_massToF2);
	    std::normal_distribution<double> gaussToF1(massToF_mean1, massToF_resolution1);
	    std::normal_distribution<double> gaussToF2(massToF_mean2, massToF_resolution2);
	    float massToF1   = gaussToF1(engine);
	    float massToF2   = gaussToF2(engine);

	    massLooseHist->Fill(massToF1);
	    massLooseHist->Fill(massToF2);

		for (int j = 0; j < nEvts_SR2Cand.size(); ++j){
			// Apply final mass requirements for each SR
			if (std::min(massToF1, massToF2)  < massToF2Cand_min[j]) continue;
			++nEvts_SR2Cand[j];
		}
		if (std::min(massToF1, massToF2)  > massToF2Cand_min[0])  ++nMassWindow;
	}
	else if (tightCandidates.size() == 1) {
		// Events with exactly one tight candiate and not two loose ones
		float mass = fabs(tightCandidates[0].m());
		// Sample the ToF mass for one full-detector candidate
		int   bin_massToF        = MToFFullDet->GetXaxis()->FindBin(mass);
		float massToF_mean       = MToFFullDet->GetBinContent(bin_massToF);
		bin_massToF        		 = MToFFullDetErr->GetXaxis()->FindBin(mass);
		float massToF_resolution = MToFFullDetErr->GetBinContent(bin_massToF);
		std::normal_distribution<double> gaussToF(massToF_mean, massToF_resolution);
		float massToF = gaussToF(engine);

	    massTightHist->Fill(massToF);

		for (int j = 0; j < nEvts_SR1Cand.size(); ++j){
			// Apply final mass requirements for each SR
			if (massToF  < massToF1Cand_min[j]) continue;
			++nEvts_SR1Cand[j];
		}
		if (massToF  > massToF1Cand_min[0]) ++nMassWindow;
	}
  // End of event loop.
  }

//  pythia.stat();

  //Save histograms to file:
  FILE* histoFile = fopen(histname.c_str(), "w");
  double xbin,ybin,binContent;
  fprintf(histoFile,"# massToF,Entries\n");
  for (int i=1; i<=massLooseHist->GetNbinsX();i++) {
	  xbin = massLooseHist->GetXaxis()->GetBinCenter(i);
	  binContent = massLooseHist->GetBinContent(i);
	  fprintf(histoFile,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fprintf(histoFile,"\n\n# massToF,Entries\n");
  for (int i=1; i<=massTightHist->GetNbinsX();i++) {
	  xbin = massTightHist->GetXaxis()->GetBinCenter(i);
	  binContent = massTightHist->GetBinContent(i);
	  fprintf(histoFile,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fprintf(histoFile,"# eta,beta,Entries\n");
  for (int i=1; i<=recHist->GetNbinsX();i++) {
	  for (int j=1; j<=recHist->GetNbinsY();j++) {
		  xbin = recHist->GetXaxis()->GetBinCenter(i);
		  ybin = recHist->GetYaxis()->GetBinCenter(j);
		  binContent = recHist->GetBinContent(i,j);
		  fprintf(histoFile,"%1.3e, %1.3e, %1.3e\n",xbin,ybin,binContent);
	  }
  }
  fprintf(histoFile,"\n\n# pT,Entries\n");
  for (int i=1; i<=pTHist->GetNbinsX();i++) {
	  xbin = pTHist->GetXaxis()->GetBinCenter(i);
	  binContent = pTHist->GetBinContent(i);
	  fprintf(histoFile,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fprintf(histoFile,"\n\n# beta,Entries\n");
  for (int i=1; i<=betaHist->GetNbinsX();i++) {
	  xbin = betaHist->GetXaxis()->GetBinCenter(i);
	  binContent = betaHist->GetBinContent(i);
	  fprintf(histoFile,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fprintf(histoFile,"\n\n# MET,Entries\n");
  for (int i=1; i<=metHist->GetNbinsX();i++) {
	  xbin = metHist->GetXaxis()->GetBinCenter(i);
	  binContent = metHist->GetBinContent(i);
	  fprintf(histoFile,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fclose(histoFile);


  fprintf(stdout,"Number of events Generated = %i\n",iEvent);
  fprintf(stdout,"Trigger eff = %1.3e\n",float(nTrigger)/float(iEvent));
  fprintf(stdout,"Decay outside eff = %1.3e\n",float(nDecay*nTrigger)/float(nHSCPs*iEvent));
  fprintf(stdout,"Candidate eff = %1.3e (average loose candidate eff = %1.3e) \n",float(nCandidate)/float(iEvent),
		  avgCandEff/float(recHist->GetEntries()));
  fprintf(stdout,"Mass window trigger eff = %1.3e\n",float(nMassWindow)/float(iEvent));

  FILE* OutputFile = fopen(outname.c_str(), "w");
  fprintf(OutputFile,"Number of events Generated = %i\n",iEvent);
  fprintf(OutputFile,"MET trigger eff = %1.3e\n",float(nTrigger)/float(iEvent));
  fprintf(OutputFile,"Decay outside eff = %1.3e\n",float(nDecay)/float(iEvent));
  fprintf(OutputFile,"Candidate eff = %1.3e (average candidate eff = %1.3e) \n",float(nCandidate)/float(iEvent),
		  avgCandEff/float(recHist->GetEntries()));
  fprintf(OutputFile,"Mass window trigger eff = %1.3e\n\n",float(nMassWindow)/float(iEvent));


  fprintf(OutputFile,"2 Candidates:\n");
  fprintf(stdout,"2 Candidates:\n");
  for (int i = 0; i < nEvts_SR2Cand.size(); ++i){
	  fprintf(OutputFile,"(mTOF>%3.0f) total Efficiency (2Cand): %1.3e +- %1.1e\n", massToF2Cand_min[i],
			float(nEvts_SR2Cand[i])/float(iEvent),sqrt(float(nEvts_SR2Cand[i]))/float(iEvent));
	  fprintf(stdout,"(mTOF>%3.0f) total Efficiency (2Cand): %1.3e +- %1.1e\n", massToF2Cand_min[i],
			float(nEvts_SR2Cand[i])/float(iEvent),sqrt(float(nEvts_SR2Cand[i]))/float(iEvent));
  }

  fprintf(OutputFile,"1 Candidate:\n");
  fprintf(stdout,"1 Candidate:\n");
  for (int i = 0; i < nEvts_SR1Cand.size(); ++i){
	  fprintf(OutputFile,"(mTOF>%3.0f) total Efficiency (1Cand): %1.3e +- %1.1e\n", massToF1Cand_min[i],
			float(nEvts_SR1Cand[i])/float(iEvent),sqrt(float(nEvts_SR1Cand[i]))/float(iEvent));
	  fprintf(stdout,"(mTOF>%3.0f) total Efficiency (1Cand): %1.3e +- %1.1e\n", massToF1Cand_min[i],
			float(nEvts_SR1Cand[i])/float(iEvent),sqrt(float(nEvts_SR1Cand[i]))/float(iEvent));
  }
  fclose(OutputFile);
    

  // Done.
  return 0;
}


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <input file>] [-o <output file>] [-d <histogram file>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <input file>:  pythia input LHE or SLHA file" << endl;
	  cout << "        -c <pythia config file>:  pythia config file [pythia8_stau.cfg]" << endl;
	  cout << "        -o <output file>:  output filename for naming the output file and histograms [<input file>.out]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]. If n < 0, it will run over all events in the LHE file" << endl;
	  cout << "        -d <histogram file>:  histogram file name. If set it will save all histograms to file [<input file>.dat]" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  float weight = 1.;
  int nevents = -1;
  double width = 0.;
  string cfgfile = "pythia8_stau.cfg";
  string outfile = "";
  string infile = "";
  string histofile = "";
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

    if ( s== "-d" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      histofile = argv[i+1];
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



    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(infile, nevents, cfgfile, outfile, histofile);

  return 0;
}
