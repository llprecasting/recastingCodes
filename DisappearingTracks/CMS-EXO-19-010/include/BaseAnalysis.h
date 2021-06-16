
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



#pragma once

#include "heputil.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <map>
#include <random>
#include <iomanip>
// Define cutflow objects and helper functions
#include "useYoda.h"
#include <cstdarg>

using namespace std;


class cut {
  public:
    double weight_left;
    int number_left;
    string name;
    bool passed;
  cut();

  cut(std::string cutname);
  ~cut();
};

class cutflow {
    private:
    
    std::map<std::string, int> cutmap;
    //double initial_weight;
    //int intial events;
    
    public:
    std::string _name;
    std::vector<cut*> cuts;
    bool passed;

    cutflow();
    ~cutflow();
    cutflow(const std::string& name);

    void addcut(const std::string& cutname);

    void applycut(const std::string& cutname,bool &condition, double &weight);

    void geteff_and_err(double &eff,double &err);
    void geteff_and_err(const std::string& cutname,double &eff,double &err);

    void print(ostream &ss);

    void print();
};



class BaseAnalysis {

private:
    //vector<YODA::AnaylsisObject*> YodaObjects;
    map<string, YODA::Histo1D*> _YodaHisto1Ds;
    map<string, YODA::Histo2D*> _YodaHisto2Ds;
    map<string, YODA::Profile1D*> _YodaProfile1Ds;
    bool _isMaster;
    
public:
    

    std::string analysisname;
    double _totalweight; 
    double _xsection;
    int _totalevents;
    HEP::Event *Event;

    std::vector<cutflow*> cutflows;
    std::map<std::string, std::vector<cutflow*>> regionmap; // maps from the name of the cut to a list of the relevant cutflows
    std::map<std::string, cutflow*> regions; // maps from the region names to the regions themselves
    double eventweight;

    
    

    virtual ~BaseAnalysis() {}; // complains if we don't define this.

    virtual void init() {}; // need to add the brackets for all virtual functions

    virtual void Execute(std::mt19937 &engine) {};

    virtual void Finalise() {};


    // To do: take this out of the init so that it only displays once per analysis ...
    //virtual void splash ();  

    void cleanup();

    void setup();

    bool ProcessEvent(HEP::Event *evnt);

    void AddRegionSelection(const std::string &region_name);

    void AddCut(const std::string& cutname,const std::string& region_name);

    void AddCut(const std::string& cutname,const std::vector<std::string> region_names);

  void ApplyCut(bool condition,const std::string& cutname);

    //void AddYodaHisto1D(const std::string& objectname, ...);
    void AddYodaHisto1D(const std::string& objectname);
    void AddYodaHisto1D(const std::string& objectname, size_t nbins, double lower, double upper);
    void AddYodaHisto1D(const std::string& objectname, const std::vector<double>& binedges);
    void FillYodaHisto1D(const std::string& objectname, double x, double fraction=1.0);
    YODA::Histo1D* GetYodaHisto1D(const std::string& objectname);

    void AddYodaHisto2D(const std::string& objectname);
    void AddYodaHisto2D(const std::string& objectname, size_t nbinsX,size_t nbinsY, double lowerX, double upperX, double lowerY, double upperY);
    void AddYodaHisto2D(const std::string& objectname, const std::vector<double>& xedges, const std::vector<double>& yedges);
    void FillYodaHisto2D(const std::string& objectname, double x, double y, double fraction=1.0);
    YODA::Histo2D* GetYodaHisto2D(const std::string& objectname);

    void AddYodaProfile1D(const std::string& objectname);
    void AddYodaProfile1D(const std::string& objectname, size_t nbins, double lower, double upper);
    void AddYodaProfile1D(const std::string& objectname, const std::vector<double>& binedges);
    void FillYodaProfile1D(const std::string& objectname, double x, double y, double fraction=1.0);
    YODA::Profile1D* GetYodaProfile1D(const std::string& objectname);

    void WriteHistos(ostream &os);

    
    void set_weight(double weight);
    void Reweight(double weight_multiplier);
    double get_weight();

    void print_cutflows();
    void set_xsection(double xs);
    

    void print_cutflows(ostream &os);
    void write_results(ostream &ss);

    bool CheckConvergence(double margin);

    //BaseAnalysis& operator = (const BaseAnalysis& B);

    //BaseAnalysis& operator += (const BaseAnalysis& B);
    void add( BaseAnalysis& B);
};
