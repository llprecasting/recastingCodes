#pragma once

#include <cassert>
#include <cstdlib>
#include <cstddef>
#include <type_traits>

#include <random>
#include <cfloat> // for DBL_MAX

//#include "time.h"
#include <chrono>
#include <sstream>
#include <string>
#include "heputil.h"
#include "include/gzstream.h"





HEP::P4 readvec(std::istringstream &ss)
{
    double x,y,z,t;
    if(! (ss >> x >> y >> z >> t))
    {
        return HEP::P4(0.0,0.0,0.0,0.0);
    }
    else
    {
        return HEP::P4(x,y,z,t);
    }


}

void ReadPileupFile(const string & minbiasfile,std::vector<HEP::PileupEvent*> *pileups, int max_events=-1)
{
    GZ::igzstream ss;
    ss.open("minbias.dat.gz");
    
    bool FoundEvent=false;
    std::string line;
    //std::stringstream line;
    int nevents=0;
    HEP::PileupEvent* new_event;
    while (std::getline(ss, line)) {
    //while (std::getline(ss, line,' ')) {
    //    std::string P_type;
     //   line >> P_type;

    //    std::cout << P_type << std::endl;
    //std::cout << line << std::endl;
    std::istringstream data(line);
    std::string Ptype;
    data >> Ptype;
    if(Ptype[0] == 'E')
    {
        
        if(nevents > 0)
        {
            
            pileups->push_back(new_event);
        }
        else
        {
            FoundEvent=true;
        }

        if((max_events>0) && (nevents > max_events)) 
        {
            ss.close();
            return;
        }
        
        new_event = new HEP::PileupEvent;
        new_event->clear();
        nevents++;

        continue;
    }
    if(!FoundEvent) continue;



    int pid,chargeType;

    if(! (data >> pid)) continue;
    if(! (data >> chargeType)) continue; 
     

    HEP::P4 mom=readvec(data);
    HEP::P4 prod=readvec(data);
    HEP::P4 dec=readvec(data);
    HEP::Particle* new_particle = new HEP::Particle(mom,pid);
    new_particle->set_prod(prod);
    new_particle->set_3Q(chargeType);

    if(dec.p2() > 0.0001)
    {
        new_particle->set_decay(dec);
    }  
    if(Ptype.compare("P") == 0)
    {
        new_event->add_particle(new_particle);

    }
    else if(Ptype.compare("CH") ==0)
    {
        new_event->add_charged_hadron(new_particle);

    }
    else if(Ptype.compare("NH")==0)
    {
        new_event->add_neutral_hadron(mom);
        delete new_particle;
    }
    else
    {
        std::cout << "Could not read line starting " << Ptype <<std::endl;
        delete new_particle;
        continue;
    }
 

    }
    ss.close();

    // add the last event!
    if(FoundEvent)
    {
        pileups->push_back(new_event);
    }
    

}

