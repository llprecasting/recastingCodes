#include "include/running.h"
#include "include/fillevent.h"


using namespace std;



void run(SettingsContainer &runSettings)
{
int nevents = runSettings.nevents;
string infile = runSettings.hepmcfile;
string outfile = runSettings.outfile;
string histofile = runSettings.histofile;
string cutflowfile = runSettings.cutflowfile;
double Rjet = runSettings.Rjet;

cout << "Finding " << nevents << endl;
 
std::srand(500); 
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 engine(seed1);

int len_filename = infile.length(); 
  
string line;

double nPileupAvg = 29;

std::vector<HEP::PileupEvent*> pileups;

if(runSettings.includePileUp)
{
  nPileupAvg=runSettings.npileup;
  ReadPileupFile(runSettings.pileupfilename,&pileups,10000);
  std::cout << "Added pileup events: " << pileups.size();
}
std::poisson_distribution<int> pileup_distribution(nPileupAvg);


// get type of file extension
std::cout << "Reading file: " << infile << std::endl;
//HepMC::IO_GenEvent* hepmcevts;
HepMC::IO_GenEvent* hepmcevts;
  // declaring character array 
char char_array[len_filename + 1]; 
// copying the contents of the 
// string to char array 
strcpy(char_array, infile.c_str()); 

std::ifstream hmcstream;
GZ::igzstream hmcgzstream;


if(infile.substr(len_filename-7,std::string::npos).find(".hepmc")!=  std::string::npos) 
{
 //hmcstream = new std::ifstream(char_array);
//hepmcevts = new HepMC::IO_GenEvent(*hmcstream);
std::cout << "Not zipped!" << std::endl;
hmcstream.open(char_array);
hepmcevts = new HepMC::IO_GenEvent(hmcstream);
}
/*
else if(infile.substr(len_filename-5,std::string::npos) == '.fifo')
{

}*/
else if(infile.substr(len_filename-10,std::string::npos).find(".hepmc.gz")!=  std::string::npos)
{
 cout << "Gzip file" << endl;

 hmcgzstream.open(char_array);
 hepmcevts = new HepMC::IO_GenEvent(hmcgzstream);
 
}
else
{
  cout << "Unrecognised file format!" << infile.substr(len_filename-9,std::string::npos) << endl;
  return;
}



 
std::cout << "Starting to read " << std::endl;




int i=0;
HepMC::GenEvent* evt = hepmcevts->read_next_event();
HepMC::Units::LengthUnit length_used=evt->length_unit();
HepMC::Units::MomentumUnit energy_used=evt->momentum_unit();
std::cout << "Read first event " << std::endl;
if(length_used != HepMC::Units::MM)
{
  cout << "HepMC file uses a length other than mm! " << length_used << " " << HepMC::Units::MM << endl;
  return;
}

if(energy_used != HepMC::Units::GEV)
{
  cout << "HepMC file uses momentum units other than GeV!" << endl;
  return;
}



double totalweight=0.0;
double survivingweight=0.0;
double evtweight;

progressbar(0.0,"");
std::string messg = "";

vector<BaseAnalysis*> Master_analyses;

for(auto analysis_name : runSettings.Analyses_To_Run)
{
  BaseAnalysis* new_analysis = mkAnalysis(analysis_name);
  new_analysis->_totalevents=0;
  new_analysis->_totalweight=0.0;
  Master_analyses.push_back(new_analysis);

}

auto starttime = std::chrono::high_resolution_clock::now();
double last_inc = 0.0;
int npileup=0;

while ( evt && (i < nevents) ) {  

  auto nowtime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timediff = std::chrono::duration_cast<std::chrono::duration<double>>(nowtime-starttime);
  double time_inc = timediff.count();
  
  std::stringstream sts;
  sts.precision(6);

  //sts << timediff.count();
  sts << (time_inc-last_inc);
  last_inc=time_inc;

  messg = sts.str();
  
  
  progressbar(((float)i)/((float) nevents),messg);
  if( (evt->weights()[0] ==0.0))
  {
  
    delete evt;
    *hepmcevts >> evt;
    i++;  
      continue; 
  }
  
  // nonzeroweight
  //evtweight=evt->weights()[0];
  //totalweight+=evtweight;
  //cutflow[0]->passed=true;
  //applycut("Initial Events",true,evtweight);
  HEP::Event HUevent;

  if(runSettings.includePileUp)
  {
    npileup=pileup_distribution(engine);
  }

  //cout << "Filling event with pileup: " << npileup << endl;
  fillevent(HUevent,evt, npileup, pileups, engine);   
  //cout << "Filled event" << endl;
  for(auto main_analysis : Master_analyses)
  {
    if(! main_analysis->ProcessEvent(&HUevent))
    {
      cout << "Event analysis failed" << endl;
      continue;
    }
    main_analysis->Execute(engine);

  }  

  HUevent.clear();

// generate next event
delete evt;
*hepmcevts >> evt;
i++;
}; // end of event loop


delete evt;
progressbar(1.0, "");
cout << endl;

for( auto Master_analysis : Master_analyses)
{
  Master_analysis->Finalise();
}


std::ofstream ResultFile,CutFlowFile,HistoFile;

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
  }
  HistoFile.close();
}

for(auto master_analysis : Master_analyses)
  {
    // delete the histograms and anything else that needs cleaning up.
    master_analysis->cleanup();
    delete master_analysis;
  }

//hmcgzstream.close();
delete hepmcevts;
// Now free up the pileup
for(auto PUEvent : pileups)
{
  delete PUEvent; //PUEvent->clear();
}



//printcutflow();
cout << "Done " << endl;
//main_analysis.~ThisAnalysis();
cout << "Returning" << endl;
  // Done.
  return;
}


void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [<yaml settings file>] [optional: hepmcfile] " << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {

  SettingsContainer runSettings;


  if(argc < 2)
  {
    help(argv[0]);
  }
 // string fname="";
 // bool foundyaml=false;
  for ( int i=1; i!=argc ; ++i )
  {
    string sss = (string) argv[i];
    if ( sss== "-h" )
    {
      help ( argv[0] );
    }

    int len_sss= sss.length();
    if(sss.substr(len_sss-6,std::string::npos).find(".yaml")!=  std::string::npos)
    {
      runSettings.ReadYAML(sss);
    }
    else
    {
      runSettings.hepmcfile=sss;
    }

  }
  //string fname=(string) argv[1];
  

  if(runSettings.Analyses_To_Run.size() ==0)
  {
    // If the user hasn't selected any analyses, add them all ...
    CreateAnalysisList(runSettings.Analyses_To_Run);
  }

  run(runSettings);  

  cout << "Finished" << endl;
  return 0;
}
