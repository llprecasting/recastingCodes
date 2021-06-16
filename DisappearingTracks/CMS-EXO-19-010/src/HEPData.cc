
#include "include/HEPData.h"
using namespace std;


void readHEPDatacsv(std::string filename, vector<vector<double>> & datas)
{
  std::string fullfilename=std::string(HEPDATAINCLUDEDIR)+filename;
  // std::cout << fullfilename << endl;
  std::ifstream ss(fullfilename);

  if(!ss.is_open())
    {
      cout << "Could not open file " << filename << endl;
      throw 1;

    }

  bool FoundStart=false;
  std::string line,cell;

  vector<double> binedges;
  vector<double> binvalues;
  datas.clear();

  
  while (std::getline(ss, line))
    {
      
      if((line[0] == '#' ) || (isalpha(line[0])))
	{
	  continue;
	}

      if(line.find(',') == std::string::npos)
	{
	  if(FoundStart)
	    {
	      // Now have come to end
	      break;
	    }
	  continue;
	}

      if(!FoundStart)
	{
	  FoundStart = true;
	}
      std::istringstream sline(line);
      vector<double> data;
      while(std::getline(sline,cell, ','))
	{
	  data.push_back(stod(cell));
	}
    // This checks for a trailing comma with no data after it.
      if (!sline && cell.empty())
	{
	  // If there was a trailing comma then add an empty element.
	  data.push_back(0.0);
	}

      datas.push_back(data);
    };

    
  // finished adding data ... Do some checks?
  
    
}


Efficiency1D::Efficiency1D(std::string filename, bool return_uncert)
{
  this->readcsv(filename,return_uncert);
};

Efficiency1D::~Efficiency1D()
{
  delete _binfn;
}

double Efficiency1D::get_at(double x)
{
  if(x < _binfn->binning.edges[0]) return 0.0;
  double rtn;
  // make safe against underflows
  try
    {
      rtn = this->_binfn->get_at(x);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;

};


void Efficiency1D::readcsv(std::string  filename, bool return_uncert)
{
  vector<vector<double>> datas;
  try
    {
      readHEPDatacsv(filename, datas);
    }
  catch (int e)
    {
      cout << "File not initialised: " << filename << endl;
      return;
    }
  // now we have our csv stored in a vector, time to process it.
  // should be of form x, xmin, xmax, y, yerror+, yerror-
  
  // first have to find the bin edges.
 
  vector<double> binedges;
  vector<double> yvals;
  for (auto dataline : datas)
    {
      if(dataline.size() < 6) continue;
       // take the lower edge ...
      binedges.push_back(dataline[1]); 
     
    }
  // For the upper range I will add infinity to avoid overflows
  binedges.push_back(DBL_MAX);
  
  // Now to define and fill the function

  HEP::BinnedFn1D<double>* _newfn = new HEP::BinnedFn1D<double>(binedges);


  if(return_uncert) // assume symmetric 
  {
    for(auto dataline : datas)
    {
      if(fabs(dataline[5])> fabs(dataline[4]))
      {
        _newfn->set_at(dataline[0],fabs(dataline[5]));
      }
      else
      {
        _newfn->set_at(dataline[0],fabs(dataline[4]));
      }
    }
  }
  else
  {
  for(auto dataline : datas)
    {
      _newfn->set_at(dataline[0],dataline[3]);
    }
  }
  this->_binfn = _newfn;

}



//////////////////////////////////////////////////////////////////////////



Efficiency2D::Efficiency2D(std::string filename, bool return_uncert)
{
  this->readcsv(filename, return_uncert);
};

Efficiency2D::~Efficiency2D()
{
  delete _binfn;
}



double Efficiency2D::get_at(double x, double y)
{
   if(x < _binfn->binning.binningX.edges[0]) return 0.0;
   if(y < _binfn->binning.binningY.edges[0]) return 0.0;
  double rtn;

  // make safe against underflows
  try
    {
      rtn = this->_binfn->get_at(x,y);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;
};



// Now for the heavy lifting

void Efficiency2D::readcsv(std::string  filename, bool return_uncert)
{
  vector<vector<double>> datas;

  readHEPDatacsv(filename, datas);
  // now we have our csv stored in a vector, time to process it.
  // should be of form x, xmin, xmax, y, ymin, ymax, zval, z+, z-
  // will ignore errors on z. 
  
  // First have to find the bin edges. I will assume square bins and that all are populated
  // If this is not the case (i.e. the binning is sparse) then we are out of luck and need
  // a tougher algorithm. 
 
  vector<double> binedgesx, binedgesy;
  vector<double> yvals;
  
  
  for (auto dataline : datas)
    {
      if(dataline.size() < 9) continue;
       // take the lower edge for both x and y ...
      double loweredgex=dataline[1];
     
     
      if((binedgesx.size() == 0 ) || (std::find(binedgesx.begin(), binedgesx.end(), loweredgex) == binedgesx.end()))
	{

	  binedgesx.push_back(loweredgex);
	}
       double loweredgey=dataline[4];
       if( (binedgesy.size() == 0 ) || (std::find(binedgesy.begin(), binedgesy.end(), loweredgey) == binedgesy.end()))
	{

	  binedgesy.push_back(loweredgey);
	}
      
     
    }

  std::sort(binedgesx.begin(), binedgesx.end());
  std::sort(binedgesy.begin(), binedgesy.end());
  // For the upper range I will add infinity to avoid overflows
  binedgesx.push_back(DBL_MAX);
  binedgesy.push_back(DBL_MAX);

  // Now to define and fill the function

  HEP::BinnedFn2D<double>* _newfn = new HEP::BinnedFn2D<double>(binedgesx, binedgesy);


  if(return_uncert) // assume symmetric 
  {
    for(auto dataline : datas)
    {
      if(fabs(dataline[8])> fabs(dataline[7]))
      {
        _newfn->set_at(dataline[0],dataline[3],fabs(dataline[8]));
      }
      else
      {
        _newfn->set_at(dataline[0],dataline[3],fabs(dataline[7]));
      }
    }
  }
  else
  {
  for(auto dataline : datas)
    {
      _newfn->set_at(dataline[0],dataline[3], dataline[6]);
    }
  }
  this->_binfn = _newfn;

}

////////////////////////////////////////////////////////////////////////////////////
// Stuff for using YODA


YodaEfficiencies::~YodaEfficiencies()
{
  for (auto hist : _1dhistos)
    {
      delete hist.second;
    }

  for (auto hist : _1dhistos_err)
    {
      delete hist.second;
    }

  for (auto hist : _2dhistos)
    {
      delete hist.second;
    }

  for (auto hist : _2dhistos_err)
    {
      delete hist.second;
    }
}


double YodaEfficiencies::get_at(std::string histoname, double x)
{

  if(_1dhistos.find(histoname) == _1dhistos.end()) return 0;
    
  if(x < _1dhistos[histoname]->binning.edges[0]) return 0.0;
  double rtn;
  // make safe against underflows
  try
    {
      rtn = _1dhistos[histoname]->get_at(x);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;

}

double YodaEfficiencies::get_at_err(std::string histoname, double x)
{

  if(_1dhistos_err.find(histoname) == _1dhistos_err.end()) return 0;
    
  if(x < _1dhistos_err[histoname]->binning.edges[0]) return 0.0;
  double rtn;
  // make safe against underflows
  try
    {
      rtn = _1dhistos_err[histoname]->get_at(x);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;

}


double YodaEfficiencies::get_at(std::string histoname, double x, double y)
{

  if(_2dhistos.find(histoname) == _2dhistos.end()) return 0;
    
  if(x < _2dhistos[histoname]->binning.binningX.edges[0]) return 0.0;
  if(y < _2dhistos[histoname]->binning.binningY.edges[0]) return 0.0;
  double rtn;
  // make safe against underflows
  try
    {
      rtn = _2dhistos[histoname]->get_at(x,y);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;

}

double YodaEfficiencies::get_at_err(std::string histoname, double x, double y)
{

  if(_2dhistos_err.find(histoname) == _2dhistos_err.end()) return 0;
    
  if(x < _2dhistos_err[histoname]->binning.binningX.edges[0]) return 0.0;
  if(y < _2dhistos_err[histoname]->binning.binningY.edges[0]) return 0.0;
  double rtn;
  // make safe against underflows
  try
    {
      rtn = _2dhistos_err[histoname]->get_at(x,y);
    }
  catch (int e)
    {
      rtn = 0.0;
    }
  return rtn;

}



HEP::BinnedFn1D<double>* YodaEfficiencies::Efficiency1D(std::string histoname)
{

  return _1dhistos[histoname];

}

HEP::BinnedFn1D<double>* YodaEfficiencies::Efficiency1D_err(std::string histoname)
{

  return _1dhistos_err[histoname];

}

HEP::BinnedFn2D<double>* YodaEfficiencies::Efficiency2D(std::string histoname)
{

  return _2dhistos[histoname];

}

HEP::BinnedFn2D<double>* YodaEfficiencies::Efficiency2D_err(std::string histoname)
{

  return _2dhistos_err[histoname];

}



void YodaEfficiencies::read(std::string filename)
{
    std::string fullfilename=std::string(HEPDATAINCLUDEDIR)+filename;
    vector<YODA::AnalysisObject*> aos = YODA::ReaderYODA::create().read(fullfilename);

  for ( YODA::AnalysisObject* ao : aos ) 
    {
      if( (ao)->type() == "Scatter2D")
	{
	  YODA::Scatter2D temp2d(*(dynamic_cast<YODA::Scatter2D*>(ao)));
	  //cout << "New: Scatter2D: " << ao->name() << endl;
	  //YODA::Histo1D  temp1d(*(dynamic_cast<YODA::Scatter2D*>(ao)));
	 

	  vector<double> binedges;
	  for(auto point : temp2d.points())
	    {
	      //double midx=point.val(0);
	      double lowx=point.val(1) - point.errMinus(1);

	      //double midy=point.val(2);
	      //double erry=std::max(point.errMinus(2),point.errPlus(2));

	      binedges.push_back(lowx);
	      
	    }

	  binedges.push_back(DBL_MAX);

	  HEP::BinnedFn1D<double>* _newfn = new HEP::BinnedFn1D<double>(binedges);
	  HEP::BinnedFn1D<double>* _newfn_err = new HEP::BinnedFn1D<double>(binedges);
	  
	  for(auto point : temp2d.points())
	    {
	      double midx=point.val(1);
	      //double lowx=point.val(0) - point.errMinus(0);

	      double midy=point.val(2);
	      double erry=std::max(point.errMinus(2),point.errPlus(2));

	      _newfn->set_at(midx,midy);
	      _newfn_err->set_at(midx,erry);
	      
	    }

	  // Now add these to the map

	  _1dhistos[ao->name()] = _newfn;
	  _1dhistos_err[ao->name()] = _newfn_err;
	  
	      
	}
      else if( (ao)->type() == "Scatter3D")
	{
	  YODA::Scatter3D temp3d(*(dynamic_cast<YODA::Scatter3D*>(ao)));
	  //cout << "New: Scatter3D: " << ao->name() << endl;

	  
	  vector<double> binedgesx,binedgesy;
	  for(auto point : temp3d.points())
	    {
	      //double midx=point.val(0);
	      double lowx=point.val(1) - point.errMinus(1);
	      double lowy=point.val(2) - point.errMinus(2);
	      //double midy=point.val(2);
	      //double erry=std::max(point.errMinus(2),point.errPlus(2));

	      if((binedgesx.size() == 0 ) || (std::find(binedgesx.begin(), binedgesx.end(), lowx) == binedgesx.end()))
		{

		  binedgesx.push_back(lowx);
		}
       
	      if( (binedgesy.size() == 0 ) || (std::find(binedgesy.begin(), binedgesy.end(), lowy) == binedgesy.end()))
		{

		  binedgesy.push_back(lowy);
		}
	      
	    }

	  std::sort(binedgesx.begin(), binedgesx.end());
	  std::sort(binedgesy.begin(), binedgesy.end());
	  
	  binedgesx.push_back(DBL_MAX);
	  binedgesy.push_back(DBL_MAX);
	  HEP::BinnedFn2D<double>* _newfn = new HEP::BinnedFn2D<double>(binedgesx,binedgesy);
	  HEP::BinnedFn2D<double>* _newfn_err = new HEP::BinnedFn2D<double>(binedgesx,binedgesy);
	  
	  for(auto point : temp3d.points())
	    {
	      double midx=point.val(1);
	      //double lowx=point.val(0) - point.errMinus(0);

	      double midy=point.val(2);

	      double val=point.val(3);
	      double err=std::max(point.errMinus(3),point.errPlus(3));

	      _newfn->set_at(midx,midy,val);
	      _newfn_err->set_at(midx,midy,err);
	      
	    }

	  // Now add these to the map

	  _2dhistos[ao->name()] = _newfn;
	  _2dhistos_err[ao->name()] = _newfn_err;

	  

	}
    }

}
