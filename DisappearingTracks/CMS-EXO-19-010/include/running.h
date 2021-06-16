



using namespace std;
#include <iostream>

#include <stdlib.h>


 


#include "heputil.h"

#include "fillPileUp.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <chrono>
#include <iomanip>

#include "include/BaseAnalysis.h"
#include "analysislist.h"


/// For setting the bar width; remove if on Windows etc
#include <sys/ioctl.h>
#include <unistd.h>
#include "yaml-cpp/yaml.h"

    
#define IF_X_RTN_CREATE_ANA_X(A)                                           \
  if (name == #A) return new A();

#define APPEND_ANALYSIS_TO_LIST(A)  \
  analysislist.push_back(#A);




BaseAnalysis* mkAnalysis(const std::string name);

void CreateAnalysisList(vector<string> &analysislist);





int shellWidth=80;

void progressbar(float progress, string message);


class SettingsContainer
{
    public:
        int ncores;
        int nevents;
        int npileup;
        string outfile;
        string cfgfile;
        string lhefile;
        string histofile;
        string cutflowfile;
        string hepmcfile;
        
        bool writehistos;
        bool includePileUp;
        bool checkconvergence;
        double convergencemargin;
        vector<std::string> Analyses_To_Run;

        vector<std::string> pythia_strings;
        string pileupfilename;
        double xs;
        double Rjet;
    

    void SetDefaults();
    SettingsContainer();
    SettingsContainer(std::string & YAMLfname);

    void ReadYAML(std::string & YAMLfname);
    template<typename T1> void SetValue(T1 & var,const string& varname, YODA_YAML::Node & config)
        {
            //string vname = (string) varname;
            YODA_YAML::Node result = config["settings"][varname];
            if(result)
            {
                var=result.as<T1>();
            }
        }
     template<typename T1> void SetValue(T1 & var,YODA_YAML::Node & result)
        {
            
            if(result)
            {
                var=result.as<T1>();
            }
        }

};



