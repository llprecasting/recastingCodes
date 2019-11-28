// Reads an input LHE or an input SLHA file and generate events using Pythia 8.
// The efficiencies for MS-Agnostic search are displayed and
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



int run(const string & infile, int nevents, const string & cfgfile, const string & outlabel, double width = 0.)
{


  std::srand(500);
  FILE* OutputFile = fopen((outlabel+".out").c_str(), "w");
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
  TH1F* MToFMSagnoErr = (TH1F*) InputFile->GetDirectory("Table 28")->Get("Hist1D_y1_e1");
  TH1F* MdedxMSagnoErr = (TH1F*) InputFile->GetDirectory("Table 27")->Get("Hist1D_y1_e1");


  // Book histograms.
  TH2F *recHist = new TH2F("Reconstruction","", 20, 0., 2.,20,0.,1.);
  TH2F *massHist = new TH2F("mRec","", 80, -2000., 4000.,80, -2000., 4000.);
  TH2F *vdecHist = new TH2F("vDec","", 100, 0., 5000., 100, 0., 10000.);
  TH1F *pTHist = new TH1F("pT","", 100, 0., 2000.);
  TH1F *betaHist = new TH1F("beta","", 100, 0., 1.);
  TH1F *metHist = new TH1F("MET","", 100, 0., 1000.);


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

  int nTrigger = 0;
  int nCandidate = 0;
  int nMassWindow = 0;
  int nDecay = 0;
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
    metHist->Fill(Etmiss);
    // Trigger decision
    if (Etmiss < 300.) {
        int   bin     = EtmissTurnOn->GetXaxis()->FindBin(Etmiss);
        float eff_Met = EtmissTurnOn->GetBinContent(bin);
        float metRandom = (std::rand()/(float)RAND_MAX);
        if (metRandom > eff_Met) continue;
//        cout << "Passed MET cut: " << metRandom << " > " << Etmiss << endl;
    }
//    cout << "Passed MET" << endl;

    ++nTrigger;

    // Events that passed the trigger
	std::vector<Particle> candidates;
	int nRhadrons = 0;
    for (int i = 0; i < event.size(); ++i) {

    	Particle particle = event[i];

		// isPrimaryRHadron identifies the initial R-Hadrons in the event record
		if (!isPrimaryRHadron(particle,event)) continue;
		++nRhadrons;
//		cout << "index = " << i << endl;
//		cout << "R-hadron = " << particle.name() << endl;
//		cout << "R-hadron charge = " << particle.charge() << endl;
//		cout << "Decay before endCal = " << decayBeforeEndHcal(particle) << endl;
//		cout << "Decay vertex = " << particle.vDec() << endl;

		// RhadCharge returns the charge of an R-hadron,
		// so we consider only R-hadrons charged after hadronisation
		if (abs(particle.charge()) != 1) continue;

		pTHist->Fill(particle.pT());
		betaHist->Fill(particle.pAbs()/particle.e());
		vdecHist->Fill(particle.vDec().pT(),fabs(particle.vDec().pz()));

		// decayBeforeEndHcal checks for a decay vertex before the end of the tile calorimeter,
		// so we only consider R-hadrons that decay after that (see measures above)
		if (decayBeforeEndHcal(particle)) continue;

		++nDecay;

//		cout << "Good candidate = " << particle.name() << " pT = " << particle.pT() << " p = " << particle.pAbs() << endl;

		// Consider only particles with minimum (transverse) momentum and maximum eta
		if (particle.pAbs() < 200.) continue;
		if (particle.pT() < 50.) continue;
		if (fabs(particle.eta()) > 1.65) continue;


		//Estimated decision
		float eta  = particle.eta();
		float beta     = particle.pAbs()/particle.e();
		int   bin_eta  = IDCaloEff->GetXaxis()->FindBin(fabs(eta));
		int   bin_beta = IDCaloEff->GetYaxis()->FindBin(beta);
		float effCand  = IDCaloEff->GetBinContent(bin_eta, bin_beta);

		avgCandEff += effCand;

		recHist->Fill(fabs(eta),beta);

		float candRandom = (std::rand()/(float)RAND_MAX);
//		cout << "beta = " << beta << " eta = " << eta << endl;
//		cout << "Calorimeter eff = " << effCand << " random = " << candRandom << endl;
		if (candRandom > effCand) continue;

		candidates.push_back(particle);
	// End of particle loop.
    }

	if (candidates.size() == 0) continue;

	++nCandidate;

	float Mass = 0.0;
	std::vector<int> passed_SR = {0,0,0,0};
	for (int i = 0; i < candidates.size(); ++i){

		Mass = fabs(candidates[i].m());
		// Events with at least one candidate
		// Sample the masses and apply the final mass window requirements
		int bin_massToF  = MToFMSagno->GetXaxis()->FindBin(Mass);
		float massToF_mean = MToFMSagno->GetBinContent(bin_massToF);
		bin_massToF  = MToFMSagnoErr->GetXaxis()->FindBin(Mass);
		float massToF_resolution = MToFMSagnoErr->GetBinContent(bin_massToF);
		std::normal_distribution<double> gaussToF(massToF_mean, massToF_resolution);
		float massToF = gaussToF(engine);


		int bin_massdEdx = MdedxMSagno->GetXaxis()->FindBin(Mass);
		float massdEdx_mean = MdedxMSagno->GetBinContent(bin_massdEdx);
		bin_massdEdx = MdedxMSagnoErr->GetXaxis()->FindBin(Mass);
		float massdEdx_resolution = MdedxMSagnoErr->GetBinContent(bin_massdEdx);
		std::normal_distribution<double> gaussdEdx(massdEdx_mean, massdEdx_resolution);
		float massdEdx   = gaussdEdx(engine);

		massHist->Fill(massdEdx,massToF);
//		cout << "Mass = " << Mass << " massToF = " << massToF << " massdEdx = " << massdEdx << endl;

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

	if (passed_SR[1] > 0) ++nMassWindow;

//    cout << "-------------- Event " << iEvent << "----------------" << endl;
  // End of event loop.
  }

//  pythia.stat();

  //Save histograms to file:
  FILE* OutputFileA = fopen((outlabel+"_massHist.dat").c_str(), "w");
  double xbin,ybin,binContent;
  fprintf(OutputFileA,"# massdEdx,massToF,Entries\n");
  for (int i=1; i<=massHist->GetNbinsX();i++) {
	  for (int j=1; j<=massHist->GetNbinsY();j++) {
		  xbin = massHist->GetXaxis()->GetBinCenter(i);
		  ybin = massHist->GetYaxis()->GetBinCenter(j);
		  binContent = massHist->GetBinContent(i,j);
		  fprintf(OutputFileA,"%1.3e, %1.3e, %1.3e\n",xbin,ybin,binContent);
	  }
  }
  fclose(OutputFileA);
  FILE* OutputFileB = fopen((outlabel+"_recHist.dat").c_str(), "w");
  fprintf(OutputFileB,"# eta,beta,Entries\n");
  for (int i=1; i<=recHist->GetNbinsX();i++) {
	  for (int j=1; j<=recHist->GetNbinsY();j++) {
		  xbin = recHist->GetXaxis()->GetBinCenter(i);
		  ybin = recHist->GetYaxis()->GetBinCenter(j);
		  binContent = recHist->GetBinContent(i,j);
		  fprintf(OutputFileB,"%1.3e, %1.3e, %1.3e\n",xbin,ybin,binContent);
	  }
  }
  fclose(OutputFileB);
  FILE* OutputFileBB = fopen((outlabel+"_vdecHist.dat").c_str(), "w");
  fprintf(OutputFileBB,"# R,|z|,Entries\n");
  for (int i=1; i<=vdecHist->GetNbinsX();i++) {
	  for (int j=1; j<=vdecHist->GetNbinsY();j++) {
		  xbin = vdecHist->GetXaxis()->GetBinCenter(i);
		  ybin = vdecHist->GetYaxis()->GetBinCenter(j);
		  binContent = vdecHist->GetBinContent(i,j);
		  fprintf(OutputFileBB,"%1.3e, %1.3e, %1.3e\n",xbin,ybin,binContent);
	  }
  }
  fclose(OutputFileBB);


  FILE* OutputFileC = fopen((outlabel+"_pTHist.dat").c_str(), "w");
  fprintf(OutputFileC,"# pT,Entries\n");
  for (int i=1; i<=pTHist->GetNbinsX();i++) {
	  xbin = pTHist->GetXaxis()->GetBinCenter(i);
	  binContent = pTHist->GetBinContent(i);
	  fprintf(OutputFileC,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fclose(OutputFileC);

  FILE* OutputFileD = fopen((outlabel+"_betaHist.dat").c_str(), "w");
  fprintf(OutputFileD,"# beta,Entries\n");
  for (int i=1; i<=betaHist->GetNbinsX();i++) {
	  xbin = betaHist->GetXaxis()->GetBinCenter(i);
	  binContent = betaHist->GetBinContent(i);
	  fprintf(OutputFileD,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fclose(OutputFileD);

  FILE* OutputFileE = fopen((outlabel+"_metHist.dat").c_str(), "w");
  fprintf(OutputFileE,"# MET,Entries\n");
  for (int i=1; i<=metHist->GetNbinsX();i++) {
	  xbin = metHist->GetXaxis()->GetBinCenter(i);
	  binContent = metHist->GetBinContent(i);
	  fprintf(OutputFileE,"%1.3e, %1.3e\n",xbin,binContent);
  }
  fclose(OutputFileE);


  fprintf(stdout,"Number of events Generated = %i\n",iEvent);
  fprintf(stdout,"MET trigger eff = %1.3e\n",float(nTrigger)/float(iEvent));
  fprintf(stdout,"Decay outside eff = %1.3e\n",float(nDecay)/float(iEvent));
  fprintf(stdout,"Candidate eff = %1.3e (average candidate eff = %1.3e) \n",float(nCandidate)/float(iEvent),
		  avgCandEff/float(recHist->GetEntries()));
  fprintf(stdout,"Mass window trigger eff = %1.3e\n",float(nMassWindow)/float(iEvent));

  fprintf(OutputFile,"Number of events Generated = %i\n",iEvent);
  fprintf(OutputFile,"MET trigger eff = %1.3e\n",float(nTrigger)/float(iEvent));
  fprintf(OutputFile,"Decay outside eff = %1.3e\n",float(nDecay)/float(iEvent));
  fprintf(OutputFile,"Candidate eff = %1.3e (average candidate eff = %1.3e) \n",float(nCandidate)/float(iEvent),
		  avgCandEff/float(recHist->GetEntries()));
  fprintf(OutputFile,"Mass window trigger eff = %1.3e\n\n",float(nMassWindow)/float(iEvent));


  for (int i = 0; i < nEvts_SR.size(); ++i){
      fprintf(OutputFile,"(mTOF>%3.0f, mdEdx > %3.0f) total Efficiency: %1.3e\n", massToF_min[i],massdEdx_min[i],
      		float(nEvts_SR[i])/float(iEvent));

      fprintf(stdout,"(mTOF>%3.0f, mdEdx > %3.0f) total Efficiency: %1.3e\n", massToF_min[i],massdEdx_min[i],
        		float(nEvts_SR[i])/float(iEvent));
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
	  cout << "        -o <output file>:  output label for naming the output file and histograms [test]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]. If n < 0, it will run over all events in the LHE file" << endl;
	  cout << "        -w <width (Gev)>: Optional width to be used [0]" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  float weight = 1.;
  int nevents = -1;
  double width = 0.;
  string cfgfile = "pythia8_gluino.cfg";
  string outlabel = "test";
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
      outlabel = argv[i+1];
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

  int r = run(infile, nevents, cfgfile, outlabel, width);

  return 0;
}
