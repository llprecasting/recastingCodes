#include "include/running.h"
#include "include/filleventP.h"

void run(SettingsContainer &runSettings)
{
int nevents = runSettings.nevents;
int numthreads =runSettings.ncores;
string cfgfile = runSettings.cfgfile;
string outfile = runSettings.outfile;
string histofile = runSettings.histofile;
string cutflowfile = runSettings.cutflowfile;
string infile = runSettings.lhefile;
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

  int trand;
  string seeds[numthreads];
  int effseeds[numthreads];
  int teffseed;


  

 //Pythia8::Pythia pythias[numthreads] { {"",false} };
 vector<Pythia8::Pythia*> pythias;
  for(int i=0; i<numthreads;i++)
      {
        pythias.push_back(new Pythia8::Pythia("",false));
	string tinfile;
	    if (numthreads ==1)
	      {
		      tinfile=infile;
	      }
	    else
	      {
		      tinfile=infile+"_"+std::to_string(i);
	      };
	        tinfile=tinfile+".lhe";
	        ifstream tempstream(tinfile);
	        if(!tempstream)
	        {
		        tinfile=tinfile+".gz";
		        ifstream tempstream(tinfile);
		        if(!tempstream)
		        {
		          cout << "Could not find input file " << tinfile << endl;
		          return;						       
		        }

	      };
		  tempstream.close();
	    trand=pythseeddistribution(seeder);
	    cout << "trand: " << trand << endl;
	    seeds[i] = std::to_string(trand);
	//teffseed= seeder()*1;
	    teffseed = EfficiencyDistribution(seeder);
	    effseeds[i] =teffseed;
	    cout << "produced effseed " << std::to_string(teffseed) << " ... " << seeder.min() << " ... " << seeder.max() << endl;


	    pythias[i]->readFile( cfgfile );
    	cout << "Using LHE file as input for thread " << i << endl;
    	pythias[i]->readString("Beams:frameType = 4");
    	pythias[i]->readString("Beams:LHEF = " + tinfile);
	    pythias[i]->readString("Random:setSeed = on");
	    cout << "seed " << seeds[i] << " for thread " << i << endl;
    	pythias[i]->readString("Random:seed = "+seeds[i]);

		
      };


  // int nMerge = pythias[0].mode("Merging:nJetMax");
  //int nMerge = pythias[0].mode("JetMatching:nJetMax");

   // Check if jet matching should be applied.
  bool doMatch   = pythias[0]->settings.flag("JetMatching:merge");

  // Check if internal merging should be applied.
  bool doMerge   = !(pythias[0]->settings.word("Merging:Process").compare("void")
    == 0);

  
  
  if (doMatch && doMerge) {
    cerr << " Jet matching and merging cannot be used simultaneously.\n" << std::endl;
    return;
  }


vector<BaseAnalysis*> Master_analyses;

for(auto analysis_name : runSettings.Analyses_To_Run)
{
  BaseAnalysis* new_analysis = mkAnalysis(analysis_name);
  new_analysis->_totalevents=0;
  new_analysis->_totalweight=0.0;
  Master_analyses.push_back(new_analysis);

}
 
 int iAbort = 0;
int iEvent=0;
auto starttime = std::chrono::high_resolution_clock::now();
 int nConvergeNext=5000;
 bool AllConverged=false;


/*
for(int i=0; i<numthreads;i++)
      {




      };

*/

#pragma omp parallel shared(iEvent,Master_analyses,AllConverged,nConvergeNext) num_threads(numthreads) //spawn the threads
{
  int threadnum = omp_get_thread_num();


//double totalweight=0.0;
//double survivingweight=0.0;



double evtweight;



Pythia8::CombineMatchingInput combined;
 if(doMatch)
   {
     combined.setHook(*pythias[threadnum]);
   }

shared_ptr<Pythia8::amcnlo_unitarised_interface> setting;

int scheme = ( pythias[threadnum]->settings.flag("Merging:doUMEPSTree")
                || pythias[threadnum]->settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( pythias[threadnum]->settings.flag("Merging:doUNLOPSTree")
                || pythias[threadnum]->settings.flag("Merging:doUNLOPSSubt")
                || pythias[threadnum]->settings.flag("Merging:doUNLOPSLoop")
                || pythias[threadnum]->settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
 //int scheme = 1; // UMEPSTree only
 //Pythia8::amcnlo_unitarised_interface* setting = NULL;
 if(doMerge)
   {
     cout << "Merging using scheme " << scheme << endl;
     //setting = new Pythia8::amcnlo_unitarised_interface(scheme);
     setting = make_shared<Pythia8::amcnlo_unitarised_interface>(scheme);
     pythias[threadnum]->setUserHooksPtr(setting);
   }

  if(doMatch)
  {
    cout << "Using matching " << endl;
  }
   
 pythias[threadnum]->init();

 // If I read in the cfg file again here it seems to screw up the CKKWL merging
 //pythias[threadnum]->readFile(cfgfile);

 // Load additional pythia options
 for(auto pythstring : runSettings.pythia_strings)
 {
   pythias[threadnum]->readString(pythstring);
 }


 //pythias[threadnum]->particleData.listAll(); 
 //pythias[threadnum]->readString("1000022:onMode = 0");
 //pythias[threadnum]->readString("1000022:tau0 = 1.0");
 //pythias[threadnum]->readString("1000022:oneChannel = 1 1.0 100 12      -13       11");

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


 
 int npileup=0;
double last_inc = 0.0;
//while((iEvent < nevents) && (!AllConverged))
//{

 while((iEvent < nevents) ) {  
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
    if(pythias[threadnum]->info.atEndOfFile()) break; 
      if (++iAbort < 50) continue;
      cout << " Event generation aborted due to too many errors!\n";
      break;
    }
  ++iEvent;
  
  Pythia8::Event& evt = pythias[threadnum]->event;

  
  evtweight=pythias[threadnum]->info.weight();
  //cout << "Raw weight " << evtweight;

  if(doMatch || doMerge)
    {
      //  double mergeweight    = pythias[threadnum].info.mergingWeight();
      // cout << "Event weight: " << evtweight << " Merge weight: " << mergeweight ;
      //	    evtweight *= mergeweight;
      //double mergeweight = pythias[threadnum]->info.mergingWeight();
	     evtweight = pythias[threadnum]->info.weightValueByIndex();
       //cout << " ... match/merged weight " << evtweight << ", " << mergeweight << endl ;
	    //cout << "  Event weight: " << evtweight << std::endl;
    }
  //cout << endl;


  if((evtweight == 0.)) continue; 
  
if(runSettings.includePileUp)
  {
    npileup=pileup_distribution(engine);
  }

  HEP::Event HUevent;
  //cout << "Callling with npileup: " << npileup << std::endl;
  filleventPYTHIA(HUevent,evt, evtweight, npileup, pileups, engine);   
  
  // ProcessEvent just sets up the cutflows and updates the weight etc
  for(auto main_analysis : main_analyses)
  {
    if(! main_analysis->ProcessEvent(&HUevent))
    {
      cout << "Event analysis failed" << endl;
      continue;
    }
    main_analysis->Execute(engine);

  }
 
  HUevent.clear();
  

 //}; // End of subrun (convergence check)

 

}; // end of event loop

// Really ought to do better with regards to the merged xsection estimation, which we should do at the beginning ...
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

cout << "Done " << endl;
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
