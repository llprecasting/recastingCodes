//Some helper functions for computing relevant quantities
#include <string>
#include <vector>
#include <vector>


using namespace Pythia8;

bool isHSCP(Particle &particle, Event &event)
//Return True (False) if the particle is (is not) a HSCP
{

	//Skip neutral particles
	if (!particle.isCharged()) return false;
	//Skip color charged particles or hadrons
	if (particle.colType() != 0) return false;
	if (particle.isHadron()) return false;
	//Skip light particles
	if (fabs(particle.m()) < 20.) return false;
	//Check if it is not an intermediate state (has itself as daughter)
	std::vector<int> daughters = particle.daughterListRecursive();
	for (int i = 0; i < daughters.size(); ++i){
		int idaughter = daughters[i];
		Particle daughterParticle = event[idaughter];
		if (daughterParticle.id() == particle.id()) return false;
	}

	return true;
}


bool isPrimaryRHadron(Particle &particle, Event &event)
//Return True (False) if the particle is (is not) a primary R-hadron
{

	if (abs(particle.status()) != 104) {return false;}
	//Check if it is not an intermediate R-hadron (has an R-hadron as daughter)
	std::vector<int> daughters = particle.daughterListRecursive();
	for (int i = 0; i < daughters.size(); ++i){
		int idaughter = daughters[i];
		Particle daughterParticle = event[idaughter];
		if (abs(daughterParticle.status()) == 104){return false;}
	}

	return true;
}

bool decayBeforeEndHcal(Particle &particle)
//Return True (False) if the particle decays before (after) the hadronic calorimeter
{
	if (particle.isFinal()){return false;}
    Vec4 decayVertex = particle.vDec();
	if (decayVertex.pT() < 3.9e3){return true;}
	if (abs(decayVertex.pz()) < 6.1e3){return true;}

	return false;
}

bool decayInsideAtlas(Particle &particle)
//Return True (False) if the particle decays inside the detector
{
	if (particle.isFinal()){return false;}
    Vec4 decayVertex = particle.vDec();
    if (decayVertex.pT() < 12e3){return true;}
    if (abs(decayVertex.pz()) < 23e3){return true;}

	return false;
}



Vec4 getMissingMomentum(Event &event)
//Returns the missing 4-vector for the event (sums up neutralinos and neutrinos)
{
    Vec4 missingETvec;
    std::vector<int> metaStableMothers;
    //First check for all final particles:
    for (int i = 0; i < event.size(); ++i) {
    	Particle ptc = event[i];
    	if (!ptc.isFinal()) continue;
    	//If particle has been produced inside or after the calorimeter, skip it and include its parent
    	if (abs(ptc.vProd().pz()) > 6100. || ptc.vProd().pT() > 2300.){
    		int mom = ptc.mother1();
    		//cout << "Adding mom " << mom << " from " << ptc.name() << " and index " << i << endl;
    		if(std::find(metaStableMothers.begin(), metaStableMothers.end(), mom) == metaStableMothers.end()) {
    			metaStableMothers.push_back(mom);
    		}
    		continue;
    	}

    	//Muons are included in MET:
    	if (ptc.idAbs() == 13){
    		//cout << "Including: " << i << " " << ptc.name() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	//Neutrinos are included in MET:
    	if (ptc.idAbs() == 12 || ptc.idAbs() == 14 || ptc.idAbs() == 16){
    		//cout << "Including: " << i << " " << ptc.name() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	if (ptc.isHadron()) continue;  //Ignore usual hadrons
    	if (ptc.idAbs() == 22) continue; //Ignore photons
    	if (ptc.isVisible() && ptc.m() < 10.) continue; //Ignore light visible particles

        //Add to missing momentum
    	//cout << "Including: " << i << " " << ptc.name() << endl;
        missingETvec += ptc.p();
    }

    //Now check all meta-stable mothers:
    for (int imom = 0; imom < metaStableMothers.size(); ++imom) {
    	Particle ptc = event[metaStableMothers[imom]];

    	//Muons are included in MET:
    	if (ptc.idAbs() == 13){
    		//cout << "Including mom: " << metaStableMothers[imom] << " " << ptc.id() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	//Neutrinos are included in MET:
    	if (ptc.idAbs() == 12 || ptc.idAbs() == 14 || ptc.idAbs() == 16){
    		//cout << "Including mom: " << metaStableMothers[imom] << " " << ptc.id() << endl;
    		missingETvec += ptc.p();
    		continue;
    	}
    	if (ptc.isHadron()) continue;  //Ignore usual hadrons
    	if (ptc.idAbs() == 22) continue; //Ignore photons
    	if (ptc.isVisible() && ptc.m() < 10.) continue; //Ignore light visible particles

        //Add to missing momentum
    	//cout << "Including mom: " << metaStableMothers[imom] << " " << ptc.id() << endl;
        missingETvec += ptc.p();
    }

    return missingETvec;
}
