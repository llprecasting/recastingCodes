#include "include/running.h"
#include "include/filleventP.h"



//const vector<std::string> Analyses_to_run = {"HSCP_ATLAS","DT_CMS"};





/*
void progressbar(float progress, string message)
{
  static int barWidth = 50;
  int pos = barWidth * progress;
  std::string lineout;

 lineout += "[";
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) lineout += "\033[1;32m=\033[0m" ;  // makes the text blue see https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
        else if (i == pos) lineout += "\033[1;32m>\033[0m";
        else lineout += " ";
    }
  lineout += "] " + to_string(int(progress * 100.0)) + " % (" + message + ")";

  for(int i =lineout.length(); i<shellWidth; i++)
  {
    lineout += " ";
  }
  
  std::cout << lineout << "\r";

  std::cout.flush();
};
*/

//void run(int nevents, const string & infile, const string & outfile, int numthreads)
void run(SettingsContainer &runSettings)
{
int nevents = runSettings.nevents;
int numthreads =runSettings.ncores;
string infile = runSettings.cfgfile;
string outfile = runSettings.outfile;
string histofile = runSettings.histofile;
string cutflowfile = runSettings.cutflowfile;
double Rjet = runSettings.Rjet;

cout << "Finding " << nevents << " across " << numthreads << " threads"<< endl;
 
struct winsize w;
ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
shellWidth=w.ws_col;

//std::srand(500); 
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 seeder(seed1);

  std::uniform_int_distribution<int> pythseeddistribution(1,900000000);
  std::uniform_int_distribution<int> EfficiencyDistribution(0,RAND_MAX);
 std::srand(EfficiencyDistribution(seeder));






 
int len_filename = infile.length(); 
  
string line;





 

/*
string file_ext="";
size_t idx = infile.rfind('.', infile.length());
   if (idx != string::npos) {
     file_ext=infile.substr(idx+1, infile.length() - i));
   }
*/

  //FILE* OutputFile = fopen(outfile.c_str(), "w");

////// PILEUP STUFF HERE!!!!   

// Data in 2016 so mu = 24.9

// see https://twiki.cern.ch/twiki/bin/view/AtlasPublic/LuminosityPublicResultsRun2#Pileup_Interactions_and_Data_AN1

double nPileupAvg = 29;

//nPileupAvg=0.0;

std::vector<HEP::PileupEvent*> pileups;

if(runSettings.includePileUp)
{
  nPileupAvg=runSettings.npileup;
  ReadPileupFile(runSettings.pileupfilename,&pileups,10000);
  std::cout << "Added pileup events: " << pileups.size();
}
std::poisson_distribution<int> pileup_distribution(nPileupAvg);

/*
ReadPileupFile("minbias.dat.gz",&pileups,10000);

std::cout << "Added pileup events: " << pileups.size();
std::poisson_distribution<int> pileup_distribution(nPileupAvg);
*/

 // pythia.readString("32:tau0 = 250.0");


  int trand;
string seeds[numthreads];
 int effseeds[numthreads];
 int teffseed;

vector<Pythia8::Pythia*> pythias;

 
 //Pythia8::Pythia pythias[numthreads] { {"",false}} ;
 //vector<Pythia8::Pythia> pythias(numthreads,Pythia8::Pythia("",false));
 
  for(int i=0; i<numthreads;i++)
      {
   //pythias.push_back(Pythia8::Pythia("",false));    
   pythias.push_back(new Pythia8::Pythia("",false));
	trand=pythseeddistribution(seeder);
	cout << "trand: " << trand << endl;
	seeds[i] = std::to_string(trand);
	//teffseed= seeder()*1;
	teffseed = EfficiencyDistribution(seeder);
	effseeds[i] =teffseed;
	cout << "produced effseed " << std::to_string(teffseed) << " ... " << seeder.min() << " ... " << seeder.max() << endl;


	pythias[i]->readFile( infile );
    	//cout << "Using LHE file as input for thread " << i << endl;
    	//pythias[i].readString("Beams:frameType = 4");
    	//pythias[i].readString("Beams:LHEF = " + tinfile);
	pythias[i]->readString("Random:setSeed = on");
	cout << "seed " << seeds[i] << " for thread " << i << endl;
	pythias[i]->readString("Random:seed = "+seeds[i]);
	
      };


  int iAbort = 0;
  //double nCuts = 0.;



vector<BaseAnalysis*> Master_analyses;

for(auto analysis_name : runSettings.Analyses_To_Run)
{
  BaseAnalysis* new_analysis = mkAnalysis(analysis_name);
  new_analysis->_totalevents=0;
  new_analysis->_totalweight=0.0;
  Master_analyses.push_back(new_analysis);

}

int iEvent=0;
auto starttime = std::chrono::high_resolution_clock::now();
 
#pragma omp parallel shared(iEvent,Master_analyses) num_threads(numthreads) //spawn the threads
{
  int threadnum = omp_get_thread_num();


double totalweight=0.0;
double survivingweight=0.0;
double evtweight;
pythias[threadnum]->init();

// Add in extra settings *after* the init
for(auto pythstring : runSettings.pythia_strings)
 {
   pythias[threadnum]->readString(pythstring);
 }

std::mt19937 engine(effseeds[threadnum]);
cout << "initialising engine with " << effseeds[threadnum] << endl;

 
progressbar(0.0,"");
std::string messg = "";


vector<BaseAnalysis*> main_analyses;

for(auto analysis_name : runSettings.Analyses_To_Run)
{
  BaseAnalysis* new_analysis = mkAnalysis(analysis_name);
  new_analysis->_totalevents=0;
  main_analyses.push_back(new_analysis);

}

/*
HSCP_ATLAS main_analysis;
main_analysis._totalevents=0;
*/

double last_inc = 0.0;

int npileup=0;

 while(iEvent < nevents) {  
//cout << evt->weights().size() << " ... " << evt->weights()[0] <<  " ... " << evt->weights()[1] <<  " ... " << evt->weights()[2] << endl;
  
  auto nowtime = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> timediff = std::chrono::duration_cast<std::chrono::duration<double>>(nowtime-starttime);
  double time_inc = timediff.count();
  
  std::stringstream sts;
  sts.precision(6);

  //sts << timediff.count();
  sts << (time_inc-last_inc);
  last_inc=time_inc;

  messg = sts.str();

  progressbar(((float)iEvent)/((float) nevents),messg);

  if (!pythias[threadnum]->next()) {
      if (++iAbort < 50) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
  ++iEvent;
  
  Pythia8::Event& evt = pythias[threadnum]->event;

  evtweight = pythias[threadnum]->info.weight();




  if(!(evtweight > 0.)) continue; 
  

  if(runSettings.includePileUp)
  {
    npileup=pileup_distribution(engine);
  }
  // otherwise have already set npileup to zero

  HEP::Event HUevent;
  //cout << "Callling with npileup: " << npileup << std::endl;
  filleventPYTHIA(HUevent,evt, evtweight, npileup, pileups, engine);   
  

  for(auto main_analysis : main_analyses)
  {
    if(! main_analysis->ProcessEvent(&HUevent))
    {
      cout << "Event analysis failed" << endl;
      continue;
    }
    main_analysis->Execute(engine);

  }
  /*
  // ProcessEvent just sets up the cutflows and updates the weight etc
  if(! main_analysis.ProcessEvent(&HUevent))
  {
    cout << "Event analysis failed" << endl;

  }
  
  //auto start_analysis = std::chrono::high_resolution_clock::now();
  //cout << "Starting analysis" << endl;
  main_analysis.Execute(engine);
  //auto end_analysis = std::chrono::high_resolution_clock::now();
  */

  //std::chrono::duration<double> analysis_diff = std::chrono::duration_cast<std::chrono::duration<double>>(end_analysis -start_analysis);
  //double timediff = difftime(endtime, starttime) * 1000.0;

  //cout.precision(10);
  //cout << "Analysis took " << analysis_diff.count() << std::endl;
  HUevent.clear();
  

};

for(int acount = 0; acount < Master_analyses.size(); acount++)
{
  if(runSettings.xs > 0.0)
  {
    Master_analyses[acount]->set_xsection(runSettings.xs);
  }
  else
  {
    Master_analyses[acount]->set_xsection(pythias[threadnum]->info.sigmaGen()*1.e12);
  }

  Master_analyses[acount]->add(*main_analyses[acount]);
  main_analyses[acount]->cleanup();
  delete main_analyses[acount];
}

delete pythias[threadnum];

}; // end of mpi pragma


for( auto Master_analysis : Master_analyses)
{
  Master_analysis->Finalise();
}


// Now free up the pileup
for(auto PUEvent : pileups)
{
  delete PUEvent; //PUEvent->clear();
}
//pythia.~Pythia();
progressbar(1.0, "");
cout << endl;

auto endtime = std::chrono::high_resolution_clock::now();

 std::chrono::duration<double> total_diff = std::chrono::duration_cast<std::chrono::duration<double>>(endtime - starttime);
cout << "Total run took " << total_diff.count() << std::endl;
  //double timediff = difftime(endtime, starttime) * 1000.0;

  //cout.precision(10);
  //cout << "Analysis took " << analysis_diff.count() << std::endl;
//main_analysis.Finalise();

//main_analysis.print_cutflows();

//master_analysis.print_cutflows();


//main_analysis.set_xsection(pythia.info.sigmaGen()*1.e12);
std::ofstream ResultFile,CutFlowFile,HistoFile;
//OutputFile.open(outfile,std::ofstream::out | std::ofstream::app);
ResultFile.open(outfile,std::ofstream::out);
CutFlowFile.open(cutflowfile,std::ofstream::out);


for(auto master_analysis : Master_analyses)
{
  master_analysis->write_results(ResultFile);

  master_analysis->print_cutflows(CutFlowFile);

}

ResultFile.close();
CutFlowFile.close();

if(runSettings.writehistos) 
{
  HistoFile.open(histofile,std::ofstream::out);

  for(auto master_analysis : Master_analyses)
  {
    master_analysis->WriteHistos(HistoFile);
    // and delete the histograms
  }
  HistoFile.close();
}

for(auto master_analysis : Master_analyses)
  {
    // delete the histograms and anything else that needs cleaning up.
    master_analysis->cleanup();
    delete master_analysis;
  }

//printcutflow();
cout << "Done " << endl;
//main_analysis.~ThisAnalysis();
//cout << "Returning" << endl;
  // Done.
 // return 0;
}



void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] <yaml input filename>" << endl;
    exit( 0 );
};

int main( int argc, const char * argv[] ) {

  SettingsContainer runSettings;
  
  runSettings.ncores = omp_get_max_threads()/2;

  if(argc < 2)
  {
    help(argv[0]);
  }
  string fname=(string) argv[1];
  runSettings.ReadYAML(fname);

  if(runSettings.Analyses_To_Run.size() ==0)
  {
    // If the user hasn't selected any analyses, add them all ...
    CreateAnalysisList(runSettings.Analyses_To_Run);
  }

  run(runSettings);
  cout << "Finished" << endl;
  return 0;


}




/*


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-i <cfg file>] [-n <number of events>] [-o <output file>] [-N <number of cores>] " << endl;
    cout << "        -i cfg file>:  configuration filename with pythia commands [test.cfg]" << endl;
	  cout << "        -o <output file>:  output file name [test.eff]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]" << endl;
    cout << "        -N <number of cores>:  Number of threads. Defaults to half of the maximum threads" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {

  SettingsContainer runSettings;
  
  runSettings.nevents = 100;
  runSettings.ncores = omp_get_max_threads()/2;
  CreateAnalysisList(runSettings.Analyses_To_Run);

  //string outfile = "test.eff";
  //string infile = "test.cfg";

  for ( int i=1; i!=argc ; ++i )
  {
    string s = argv[i];
    if ( s== "-h" )
    {
      help ( argv[0] );
    }

    if ( s== "-n" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      runSettings.nevents = atoi(argv[i+1]);
      i++;
      continue;
    }

    if ( s== "-i" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      runSettings.cfgfile = argv[i+1];
      i++;
      continue;
    }
    
    if ( s== "-o" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      runSettings.outfile = argv[i+1];
      i++;
      continue;
    }

    if ( s== "-N" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      runSettings.ncores = atoi(argv[i+1]);
      i++;
      continue;
    }

    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  run(runSettings);
  cout << "Finished" << endl;
  return 0;
}
*/