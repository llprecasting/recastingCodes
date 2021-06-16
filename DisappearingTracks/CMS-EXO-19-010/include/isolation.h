#pragma once

#include "heputil.h"
#include "math.h"
double d0Calc(const HEP::P4 &xv, const HEP::P4 &p);

double dzCalc(const HEP::P4 &xv, const HEP::P4 &p);

double dz_sintheta(const HEP::P4 &xv, const HEP::P4 &p);

//template<typename T1, typename T2> std::vector<T1*> overlapRemoval

//isolation 

double sumpTisolation(const HEP::P4 &p4, const HEP::Event* evt, const double &DR);

double sumETIsolation(const HEP::P4 &p4, const HEP::Event* evt, const double &DR);

double totalEventE(const HEP::Event* evt);
double medianEventE(const HEP::Event* evt);

double sumCaloE(const HEP::P4 &p4, const HEP::Event* evt, const double &DR);
double coneAreaApprox(const HEP::P4 &p4, const double &DR);

bool IsLooseTightLepton(const  HEP::Particle* Lep, const HEP::Event* event,
		   const double DeltaRp, const double DeltaRE, const double pratio, const double Eratio);

bool IsHighPTCaloOnly(const HEP::Particle* Lep, const HEP::Event* event, const double DeltaRE);



template<typename T1> void filterPhaseSpace(std::vector<T1*> &vec_t1, const double &pTmin, const double &absEtaMax)
{
	// filters the particles/jets in a vector by pT and eta.

  
  auto it = vec_t1.begin();
  while(it != vec_t1.end())
  {
	  //T1* t1 = vec_t1[it];
	  //if( (t1->pT() > pTmin) || (t1->abseta() > absEtaMax))
    if( ((*it)->pT() < pTmin ) || ((*it)->abseta() > absEtaMax))
	  {
       it = vec_t1.erase(it);
        }
    else
      {
	it++;
      }
    
  }
  
}

template<typename T1> void delete_and_filterPhaseSpace(std::vector<T1*> &vec_t1, const double &pTmin, const double &absEtaMax)
{
	// filters the particles/jets in a vector by pT and eta.

  std::vector<T1*> out_vector;
  for(auto t1 : vec_t1)
  {
    if( (t1->pT() < pTmin ) || (t1->abseta() > absEtaMax))
	  {
       delete t1;
        }
    else
      {
        out_vector.push_back(t1);
      }

  }
  vec_t1=out_vector;
    
}


// Overlap Removal
//template<typename T1, typename T2> std::vector<T1*> Removal(std::vector<T1*> &v1, std::vector<T2*> &v2,const double &drmin);

// Overlap Removal for electrons with Delta R < min(0.4, 0.04 + 10 GeV/pT) around Jet
//template<typename T1, typename T2> std::vector<T1*> RemovalJE(std::vector<T1*> &v1, std::vector<T2*> &v2);


// Overlap Removal

template<typename T1, typename T2> std::vector<T1*>
  FullRemoval(std::vector<T1*> &v1, std::vector<T2*> &v2,
  const double &drmin)
{
  // Determining with objects should be removed, and free memory
  
  std::vector<bool> mask(v1.size(),false);

  for (unsigned int j=0;j<v1.size();j++)
  {
    for (unsigned int i=0; !mask[j] && i<v2.size();i++)
        {
            if (v2[i]->mom().deltaR_eta(v1[j]->mom()) < drmin) 
	            {
	                mask[j]=true;
	                //break;
	            }
        }
    };
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
  {
    if (!mask[i]) { 
      cleaned_v1.push_back(v1[i]);
      }
      else
      {
        delete v1[i];
      }

  }

  return cleaned_v1;
};




template<typename T1, typename T2> std::vector<T1*>
  Removal(std::vector<T1*> &v1, std::vector<T2*> &v2,
  const double &drmin)
{
  // Determining with objects should be removed
  
  std::vector<bool> mask(v1.size(),false);

  for (unsigned int j=0;j<v1.size();j++)
  {
    for (unsigned int i=0; !mask[j] && i<v2.size();i++)
        {
            if (v2[i]->mom().deltaR_eta(v1[j]->mom()) < drmin) 
	            {
	                mask[j]=true;
	                //break;
	            }
        }
    };
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
};


template<typename T1> std::vector<T1*> SelfRemoval(std::vector<T1*> &v1, const double &drmin)
{
  // Determining with objects should be removed -- keep the higher pT version. This only applies for electrons anyway
  
  std::vector<bool> mask(v1.size(),false);

      double tdr;
      for (unsigned int j=0;j<v1.size();j++){
	      for (unsigned int i=j+1;i<v1.size();i++) {

	        tdr=v1[i]->mom().deltaR_eta(v1[j]->mom());
	        if ((tdr < drmin)) 
		      {
           if(v1[i]->mom().pT() < v1[j]->mom().pT())
            {
		          mask[i]=true;
             }
            else
            {
              mask[j]=true;
            }
		      }

	      }
      };
   

  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
};



// Overlap Removal for electrons with Delta R < min(0.4, 0.04 + 10 GeV/pT) around Jet
template<typename T1, typename T2> std::vector<T1*>
  RemovalJE(std::vector<T1*> &v1, std::vector<T2*> &v2)
{
  // Determining with objects should be removed
  double EPT,drmin;
  std::vector<bool> mask(v1.size(),false);
  for (unsigned int j=0;j<v1.size();j++)
    {
      EPT=v1[j]->pT();
      drmin=0.04+10.0/EPT;
      if(drmin > 0.4) drmin=0.4;
    
      for (unsigned int i=0;!mask[j] && i<v2.size();i++)
	      if (v2[i]->mom().deltaR_eta(v1[j]->mom()) < drmin)
	    {
        
	        mask[j]=true;
	        //break;
	    }
    }
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
}
