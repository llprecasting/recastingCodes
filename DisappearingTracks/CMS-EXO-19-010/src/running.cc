// Routines for handling initialisation and the progress bar.
// Stick all settings in here



#include "include/running.h"


// So that we can define the analyses to include dynamically and choose which to include at runtime
    BaseAnalysis* mkAnalysis(const std::string name)
    {
      MAP_ANALYSES(IF_X_RTN_CREATE_ANA_X)

      throw std::runtime_error("The analysis " + name + " is not a known analysis.");
      return nullptr;
    }

void CreateAnalysisList(vector<string> &analysislist)
{
    // This function adds *all* available analyses to the run 
    analysislist.clear();
    MAP_ANALYSES(APPEND_ANALYSIS_TO_LIST)
}



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


/////////////////////////////////////////////////
// Settings stuff




SettingsContainer::SettingsContainer()
{
    this->SetDefaults();
}

void SettingsContainer::SetDefaults()
{
    ncores=1;
    nevents=1000000; // million events should be enough as a default for most files
    xs=-1.0;
    Rjet=0.4;
    lhefile="";
    outfile="test.slha";
    cfgfile="test.cfg";
    histofile="test.yoda";
    cutflowfile="test.eff";
    hepmcfile="test.hepmc";
    writehistos=true;
    includePileUp=false;
    npileup=29;
    pileupfilename="minbias.dat.gz";
    Analyses_To_Run = {};
    checkconvergence=false;
    convergencemargin=0.1;
}





void SettingsContainer::ReadYAML(std::string & YAMLfname)
{
    this->SetDefaults();
    YODA_YAML::Node config = YODA_YAML::LoadFile(YAMLfname);
    for (const auto& it : config["analyses"])
   {
       Analyses_To_Run.push_back(it.as<std::string>());
       cout << "Adding analysis " << it.as<std::string>() << endl;  
   }

  this->pythia_strings.clear();
   for(const auto& it: config["Pythia Extras"])
   {
     pythia_strings.push_back(it.as<std::string>());
   }

    this->SetValue(this->ncores,"cores",config);
    this->SetValue(this->cfgfile,"Config file",config);
    this->SetValue(this->lhefile,"LHE file",config);
    this->SetValue(this->hepmcfile,"HEPMC file",config);
    this->SetValue(this->xs,"Cross-section",config);
    this->SetValue(this->Rjet,"Rjet",config);
    this->SetValue(this->nevents,"nevents",config);
    this->SetValue(this->includePileUp,"Include Pileup",config);
    this->SetValue(this->npileup,"Average Pileup",config);
    this->SetValue(this->checkconvergence,"Check Convergence",config);
    this->SetValue(this->convergencemargin,"Convergence Margin",config);

    this->SetValue(this->pileupfilename,"Pileup Filename",config);
    this->SetValue(this->writehistos,"Write Histograms",config);
    this->SetValue(this->outfile,"Efficiency Filename",config);
    this->SetValue(this->cutflowfile,"Cutflow Filename",config);
    this->SetValue(this->histofile,"Histogram Filename",config);

}
