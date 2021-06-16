/*

  Routines for reading HEPData csv and YODA files into efficiency functions

 */



#pragma once

#include <cmath>
#include <ctype.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>  // for find
#include <cfloat> // for DBL_MAX
#include "ExHEPUtils/BinnedFn.h"
#include "useYoda.h"


#ifndef HEPDATAINCLUDEDIR
#define HEPDATAINCLUDEDIR ""
#endif


using namespace std;


void readHEPDatacsv(std::string filename, vector<vector<double>> & datas);

  

class Efficiency1D
{
private:
  HEP::BinnedFn1D<double>*  _binfn;
  
public:

  
  
  Efficiency1D() {};
  ~Efficiency1D();

  Efficiency1D(std::string filename, bool return_uncert=false);
  
  void readcsv(std::string filename, bool return_uncert=false);

  double get_at(double x);


};



class Efficiency2D
{
private:
  HEP::BinnedFn2D<double>*  _binfn;
  
public:

  Efficiency2D() {};
  ~Efficiency2D();
        
  Efficiency2D(std::string filename, bool return_uncert=false);
  
  void readcsv(std::string filename, bool return_uncert=false);

  double get_at(double x, double y);


};



class YodaEfficiencies
{
    private:
        //map<std::string, YODA:Histo1D*> _1dhistos;
  //vector<HEP::BinnedFn1D<double>*> _1dhistos;
  //vector<HEP::BinnedFn1D<double>*> _1dhistos_err;
        map<std::string, HEP::BinnedFn1D<double>*> _1dhistos;
	map<std::string, HEP::BinnedFn1D<double>*> _1dhistos_err;

	map<std::string, HEP::BinnedFn2D<double>*> _2dhistos;
	map<std::string, HEP::BinnedFn2D<double>*> _2dhistos_err;

	
    public:
        YodaEfficiencies() {};
        ~YodaEfficiencies();
        
        void read(std::string filename);

        double get_at(std::string histoname, double x);
	double get_at_err(std::string histoname, double x);
	double get_at(std::string histoname, double x, double y);
	double get_at_err(std::string histoname, double x, double y);

	HEP::BinnedFn1D<double>* Efficiency1D(std::string histoname);
	HEP::BinnedFn1D<double>* Efficiency1D_err(std::string histoname);

	HEP::BinnedFn2D<double>* Efficiency2D(std::string histoname);
	HEP::BinnedFn2D<double>* Efficiency2D_err(std::string histoname);

};

