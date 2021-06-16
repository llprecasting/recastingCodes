/*

// My implementation of an analysis container. Written from scratch.
// Functions similarly to MadAnalysis analyses, so we have to define
// void init(), void Execute(std::mt19937 &engine), void Finalise() for a 
// new analysis *but* I also require a constructor to be defined (with some standard stuff).
// 
//
// Signal regions defined via:
// AddRegionSelection("region name");
// or 
// AddRegionSelection(<vector of names>);

// Cuts are defined via
// AddCut("cut name",<signal region>);
// or 
// AddCut("cut name",<vector of region names>);

// Cuts are applied by
// ApplyCut(bool condition,const std::string& cutname);

// Histogramming is done here via Yoda, so AddYodaHist1D in the init, FillYodaHist1D in the analysis
// etc.


*/



#include "include/BaseAnalysis.h"


cut::cut()
  {
    weight_left=0.0;
    number_left=0;
    name="";
    
  };

 cut::cut(std::string cutname)
  {
    weight_left=0.0;
    number_left=0;
    name=cutname;
    
  }



cut::~cut() { };

cutflow::~cutflow() { for(auto cut : cuts) { delete cut; }; };

cutflow::cutflow() {
        _name = "" ;
         passed=true;
         cut* newcut = new cut("Initial Events");
         cuts.push_back(newcut);
    };
cutflow::cutflow(const std::string& name)
    {
        _name=name;
        passed=true;
        cut* newcut = new cut("Initial Events");
         cuts.push_back(newcut);
    };

void cutflow::addcut(const std::string& cutname)
    {
        cut* newcut = new cut(cutname);
        cutmap[cutname]=cuts.size(); // do this before pushing back ...
        cuts.push_back(newcut);
        
    };

void cutflow::applycut(const std::string& cutname,bool &condition, double &weight)
    {
        if(!this->passed) // "this" not strictly necessary. 
        {
            return;
        }
        int cutindex=cutmap[cutname];
        if (condition)
        {
            cuts[cutindex]->weight_left+=weight;
            cuts[cutindex]->number_left+=1;
        }
        else
        {
            this->passed=false;
        }

    };

void cutflow::print(ostream &ss)
    {
        double initialweight=cuts[0]->weight_left;
        double eff=0.0;
        double Nevents=(double) cuts[0]->number_left;
        
        ss.precision(6);
        
        for(auto cut: cuts)
        {
            eff=(cut->weight_left)/initialweight;
            double reluncert=0.0;
            double uncert=0.0;
            if((Nevents > 0.0) && (eff > 0.0) && (eff <= 1.0))
            {
                uncert=sqrt(eff*(1.0-eff)/Nevents );
                reluncert=uncert/eff;
            }
            
            ss << cut->name << " : " << cut->number_left << ", " << cut->weight_left << " -> efficiency: " << eff << " +/- " << uncert << endl;

        }

    };

void cutflow::print()
    {
        print(std::cout);
    }

void cutflow::geteff_and_err(double &eff,double &err)
{
    geteff_and_err(this->cuts.back()->name,eff,err);
}
void cutflow::geteff_and_err(const std::string& cutname,double &eff,double &err)
{
    int cutindex=cutmap[cutname];
    double initialweight=cuts[0]->weight_left;
    int goodevents=cuts[0]->number_left;
    double finalweight=cuts[cutindex]->weight_left;
    eff=0.0;    
    if(initialweight>0.0) eff=finalweight/initialweight;
    double uncert=0.0;
    err=1.0; // default value should be order 1 uncertainty, basically if the efficiency is zero.
    if((goodevents > 0) && (eff > 0.0) && (eff < 1.0))
        {
            uncert=sqrt(eff*(1.0-eff)/goodevents);
            err=uncert/eff;
        } 
}

bool BaseAnalysis::CheckConvergence(double margin)
{
    //bool isConverged=true;
    double eff, err;
    for(auto region : regions)
    {
        region.second->geteff_and_err(eff,err);
        if(err > margin)
        {
            return false;
        }
    }
    return true;
}

void BaseAnalysis::cleanup()
{
    
    for(auto cf : cutflows)
    {
       // cf->~cutflow();
        delete cf; // calls destructor, then deallocates
    }
    
    

   for(auto histo : _YodaHisto1Ds)
   {
       delete histo.second;
   }

   for(auto histo : _YodaHisto2Ds)
   {
       delete histo.second;
   }
   for(auto histo : _YodaProfile1Ds)
   {
       delete histo.second;
   }

}

void BaseAnalysis::setup()
{
    _totalevents = 0;
    _xsection = 0.0;
    _totalweight = 0.0;
    eventweight=0.0;
    cutflows.clear();
    regions.clear();
    regionmap.clear();
    analysisname="";
    //this->_isMaster=isMaster;
    //_isMaster=true; // by default make it a "master" analysis, can optionally set as a "child"
}


void BaseAnalysis::set_xsection(double xs)
{
    _xsection=xs;
}

void BaseAnalysis::write_results(ostream &ss)
{
    ss << "BLOCK PROCESS " << analysisname << std::endl;
    //ss << "XS  " << pythia.info.sigmaGen()*1.e12 << "   # (fb) " << std::endl;
    ss << "XS   " << this->_xsection <<  "   # (fb) " << std::endl;
    ss << "Events  " << _totalevents << std::endl;

    ss << "BLOCK EFFICIENCIES " << analysisname << " # name, eff, rel. uncert." << std::endl;

    for(auto region : regions )
    {
         cutflow* thiscutflow=region.second;
         double initialweight=thiscutflow->cuts[0]->weight_left;
         int goodevents=thiscutflow->cuts[0]->number_left;
         double finalweight=thiscutflow->cuts.back()->weight_left;
         double eff=0.0;
         if(initialweight>0.0) eff=finalweight/initialweight;

    
            double uncert=0.0;
            double reluncert=0.0;
        if((goodevents > 0) && (eff > 0.0) && (eff < 1.0))
        {
        uncert=sqrt(eff*(1.0-eff)/goodevents);
        reluncert=uncert/eff;
        }

        ss << thiscutflow->_name << "   " << scientific <<  eff << "   " << reluncert << std::endl;
    }

}



bool  BaseAnalysis::ProcessEvent(HEP::Event *evnt){

        eventweight=evnt->weights()[0];
       
        if(eventweight == 0.0)
        {
            return false;
        }
        
        _totalweight+=eventweight;
        _totalevents++;
        Event = evnt;

        for( auto region : regions)
        {
            region.second->passed=true; // second refers to second term in iterator
            region.second->cuts[0]->weight_left+=eventweight;
            region.second->cuts[0]->number_left+=1;
        }

        return true;

    };

void BaseAnalysis::AddRegionSelection(const std::string &region_name) {
        cutflow* newreg=new cutflow(region_name);
        this->cutflows.push_back(newreg);
        regions[region_name] = newreg;
    };


void BaseAnalysis::AddCut(const std::string& cutname,const std::vector<std::string> region_names) {
    
    for(auto name: region_names)
    {
        
        this->AddCut(cutname,name);

    }
};


void BaseAnalysis::AddCut(const std::string& cutname,const std::string& region_name)
    {
        regions[region_name]->addcut(cutname);
        regionmap[cutname].push_back(regions[region_name]); // pointer to the cutflow
    }

    

void BaseAnalysis::ApplyCut(bool condition, const std::string& cutname)
    {
        vector<cutflow*> region_names=regionmap[cutname];
        for ( cutflow* region : region_names)
        {
            //std::cout << "Applying cut with weight " << this->eventweight << std::endl;
            region->applycut(cutname,condition,this->eventweight);
        }

    }



void BaseAnalysis::print_cutflows()
{
    BaseAnalysis::print_cutflows(std::cout);

}

void BaseAnalysis::print_cutflows(ostream &os)
{

    
    for( auto region : regions)
    {
        os << "============================= " << region.first << " =============================" << endl;
        region.second->print(os);
        
    }
    os << "==========================================================" << endl;
}

void BaseAnalysis::set_weight(double weight)
{
    this->eventweight = weight;

    // NO! Don't want to mess with the underlying event. 
    // We don't use it for the cutflows and
    // it may be used by other analyses
    //Event->set_weight(weight);


}

void BaseAnalysis::Reweight(double weight_multiplier)
{
    this->eventweight *= weight_multiplier;

    // NO! Don't want to mess with the underlying event. 
    // We don't use it for the cutflows and
    // it may be used by other analyses
    //Event->set_weight(this->eventweight);


}

    
double BaseAnalysis::get_weight()
{
    return this->eventweight;
}


    void BaseAnalysis::add ( BaseAnalysis& B)
    {   
        this->_totalweight+=B._totalweight;
        this-> _totalevents+=B._totalevents;
        for(auto region : this->regions )
        {
         cutflow* thiscutflow=region.second;
         
         
         cutflow* thatcutflow=B.regions[region.first];
         
            for(int cc=0; cc< thatcutflow->cuts.size(); cc++)
            {
                thiscutflow->cuts[cc]->weight_left+=thatcutflow->cuts[cc]->weight_left;
                thiscutflow->cuts[cc]->number_left+=thatcutflow->cuts[cc]->number_left;
            }
            
        
        }

        for(auto histomap : _YodaHisto1Ds)
        {
            *(histomap.second) += *(B.GetYodaHisto1D(histomap.first));

        }

        for(auto histomap : _YodaHisto2Ds)
        {
            *(histomap.second) += *(B.GetYodaHisto2D(histomap.first));

        }

        for(auto histomap : _YodaProfile1Ds)
        {
            *(histomap.second) += *(B.GetYodaProfile1D(histomap.first));

        }
       

    }




//////////////////////////////////////////////////////////////////////
// YODA Stuff
//

void BaseAnalysis::AddYodaHisto1D(const std::string& objectname)
{
    // Adds an empty Histo1D. The user can then manipulate it at will using the 
    // GetYodaHisto1D method ... or use my handily written fill function.

    //static string nullpath="";
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    string path=this->analysisname+"/"+objectname;

    YODA::Histo1D* tHist = new YODA::Histo1D(path, objectname);
    _YodaHisto1Ds[objectname] = tHist;

}
void BaseAnalysis::AddYodaHisto1D(const std::string& objectname, size_t nbins, double lower, double upper)
{
    
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    string path=this->analysisname+"/"+objectname;
    YODA::Histo1D* tHist = new YODA::Histo1D(nbins, lower, upper, path, objectname);
    _YodaHisto1Ds[objectname] = tHist;

}

void BaseAnalysis::AddYodaHisto1D(const std::string& objectname, const std::vector<double>& binedges)
{
    
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    string path=this->analysisname+"/"+objectname;
    YODA::Histo1D* tHist = new YODA::Histo1D(binedges, path, objectname);
    _YodaHisto1Ds[objectname] = tHist;

}

void BaseAnalysis::FillYodaHisto1D(const std::string& objectname, double x, double fraction)
{
    _YodaHisto1Ds[objectname]->fill(x,eventweight,fraction);
}

YODA::Histo1D* BaseAnalysis::GetYodaHisto1D(const std::string& objectname)
{
    return _YodaHisto1Ds[objectname];
}

///////////////////////////// profile1d
void BaseAnalysis::AddYodaProfile1D(const std::string& objectname)
{
   
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Profile1D* tHist = new YODA::Profile1D(objectname, objectname);
    _YodaProfile1Ds[objectname] = tHist;

}
void BaseAnalysis::AddYodaProfile1D(const std::string& objectname, size_t nbins, double lower, double upper)
{
    
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Profile1D* tHist = new YODA::Profile1D(nbins, lower, upper, objectname, objectname);
    _YodaProfile1Ds[objectname] = tHist;

}

void BaseAnalysis::AddYodaProfile1D(const std::string& objectname, const std::vector<double>& binedges)
{
    
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Profile1D* tHist = new YODA::Profile1D(binedges, objectname, objectname);
    _YodaProfile1Ds[objectname] = tHist;

}

void BaseAnalysis::FillYodaProfile1D(const std::string& objectname, double x, double y, double fraction)
{
    _YodaProfile1Ds[objectname]->fill(x,y,eventweight,fraction);
}

YODA::Profile1D* BaseAnalysis::GetYodaProfile1D(const std::string& objectname)
{
    return _YodaProfile1Ds[objectname];
}

///////////////////////////// end profile1d


///////////////////////////// histo2d
void BaseAnalysis::AddYodaHisto2D(const std::string& objectname)
{
    //static string nullpath="";
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Histo2D* tHist = new YODA::Histo2D(objectname, objectname);
    _YodaHisto2Ds[objectname] = tHist;

}
void BaseAnalysis::AddYodaHisto2D(const std::string& objectname, size_t nbinsX,size_t nbinsY, double lowerX, double upperX, double lowerY, double upperY)
{
    //static string nullpath="";
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Histo2D* tHist = new YODA::Histo2D(nbinsX,nbinsY,lowerX,upperX,lowerY,upperY, objectname, objectname);
    _YodaHisto2Ds[objectname] = tHist;

}

void BaseAnalysis::AddYodaHisto2D(const std::string& objectname, const std::vector<double>& xedges, const std::vector<double>& yedges)
{
    //static string nullpath="";
    // First copy of "objectname" here is the "path". I just set this to be the same as the object name, so that when
    // plotting it is produced in the same directory
    YODA::Histo2D* tHist = new YODA::Histo2D(xedges, yedges, objectname, objectname);
    _YodaHisto2Ds[objectname] = tHist;

}

void BaseAnalysis::FillYodaHisto2D(const std::string& objectname, double x, double y, double fraction)
{
    _YodaHisto2Ds[objectname]->fill(x,y,eventweight,fraction);
}

YODA::Histo2D* BaseAnalysis::GetYodaHisto2D(const std::string& objectname)
{
    return _YodaHisto2Ds[objectname];
}

///////////////////////////// end Histo2d




void BaseAnalysis::WriteHistos(ostream &os)
{
    for(auto histomap : _YodaHisto1Ds)
        {
            YODA::WriterYODA::write(os, *(histomap.second));

        }

    for(auto histomap : _YodaHisto2Ds)
        {
            YODA::WriterYODA::write(os, *(histomap.second));

        }

    for(auto histomap : _YodaProfile1Ds)
        {
            YODA::WriterYODA::write(os, *(histomap.second));

        }

}
