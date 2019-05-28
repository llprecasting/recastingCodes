/* Toy detector simulation for use with Pythia8
   Author: Nishita Desai (with input from Jad Marrouche)

   Jet energy resolution: Use 20% for jets at 50 GeV, falling linearly
   to 10% at 100 GeV, then flat 10%. 

   Jet energy scale: For jets with |eta| > 2, use 3% flat, for jets
     with |eta| < 2, use 1% flat (I am assuming the jets are above
     20-30 GeV by which point this is probably quite accurate)


   Electron resolution: use 2% at 10 GeV, falling linearly to 1% at
     100GeV, and then 1% flat.  Electron scale is effectively 0 so we
     can forget it.

   Muon resolution: between eta -2 and 0, linearly falling from 4% to
    1.5%. Symmetric; Muon Scale: small enough to ignore. 

*/

#include "ToyDetector-ATLAS.h"

namespace Pythia8{

  
  // Fastjet analysis - select algorithm and parameters

  double   ToyDetector::Rparam = 0.4;
  fastjet::Strategy ToyDetector::strategy = fastjet::Best;
  fastjet::RecombinationScheme ToyDetector::recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition * ToyDetector::jetDef  = new
    fastjet::JetDefinition(fastjet::kt_algorithm,
			   ToyDetector::Rparam,
			   ToyDetector::recombScheme,
			   ToyDetector::strategy);

bool ToyDetector::getObjects(const Event& event){

    if(event.size() < 1) {
      cout << "Invalid event" << endl;
      return false;
    }

    // Initialise all vectors to zero

    jets.clear();
    electrons.clear();
    muons.clear();
    leptons.clear();

    MET = Vec4(0.0, 0.0, 0.0, 0.0);

    // Reset Fastjet input
    fjInputs.resize(0);

    // Find isolated leptons, jets and smear momenta

    for (int i = 0; i < event.size(); i++){

      if (!event[i].isFinal()) continue;

      // No neutrinos
      if (event[i].idAbs() == 12 || event[i].idAbs() == 14 ||
          event[i].idAbs() == 16)     continue;

      // No uncharged stable BSM particles
      if (abs(event[i].id()) > 100000 && !event[i].isCharged()) continue;

      // Only |eta| < 4.5
      if (fabs(event[i].eta()) > 4.5) continue;

      // Production vertex within detector bounds (using cylindrical box r = 10m, z = 20m)
      double decayR = sqrt(pow2(event[i].xProd()) + pow2(event[i].yProd()));
      double decayZ = abs(event[i].zProd());

      if (decayR > 1100.0 || decayZ > 2800.0) continue;


      // Check if it is an isolated lepton
      if(abs(event[i].id()) == 11 || abs(event[i].id()) == 13) {
      
	if(abs(event[i].eta()) > 2.5 || event[i].pT() < 10.0) continue;

      
	double res = 0.0;
	double energy = event[i].e();
	int lid = abs(event[i].id());
      
	Vec4 v1 = event[i].p();
	double ptsum = 0.0;
	// Is it isolated?
	for(int i2 = 0; i2 < event.size(); i2++){
	  if(!event[i2].isFinal() || i2 == i || !event[i].isCharged()) continue;
	  Vec4 v2 = event[i2].p();
	  if(REtaPhi(v1, v2) < 0.2) ptsum += event[i2].pT();
	}
	
	// TODO: implement different isolation criteria
	
	if(lid == 11 && ptsum/event[i].pT() > 0.1) continue;
	if(lid == 13 && ptsum > 1.8) continue;
	
	// Smearing
	if(abs(event[i].id()) == 11) {
	  if(energy < 100.0)
	    res = (max(-0.00011 * energy, -0.01) + 0.021) * energy;
	  else
	    res = 0.01 * energy;
	}
	else {
	  double eta = event[i].eta();
	  if(abs(eta) < 2.0)
	    res = (min(0.0175 * energy, 0.025) + 0.015) * energy;
	  else
	    continue;
	}
      
	double enew = energy + res * rndm.gauss();
	Vec4 pnew = enew/energy * event[i].p();
      
	// Case isolated electron
	if(abs(event[i].id()) == 11) {
	  bool placed = false;
	  if(electrons.size() > 0){
	    for(unsigned long l = 0; l < electrons.size(); l++){
	      if(electrons.at(l).pT() < pnew.pT()) {
		electrons.insert(electrons.begin() + l, pnew);
		placed = true;
		break;
	      }
	    }
	  }
	  else
	    electrons.push_back(pnew);
	}
	else { 	// Case isolated muon
	  
	  bool placed = false;
	  if(muons.size() > 0){
	    for(unsigned long l = 0; l < muons.size(); l++){
	      if(muons.at(l).pT() < pnew.pT()) {
		muons.insert(muons.begin() + l, pnew);
		placed = true;
		break;
	      }
	    }
	  }
	  
	}
      
	// Default fill in leptons
	
	bool placed = false;
	if(leptons.size() > 0){
	  for(unsigned long l = 0; l < leptons.size(); l++){
	    if(leptons.at(l).pT() < pnew.pT()) {
	      leptons.insert(leptons.begin() + l, pnew);
	      placed = true;
	      break;
	    }
	  }
	}
      
	MET -= pnew;
	// End isolated leptons
      } else {

	// Store as input to Fastjet
	fjInputs.push_back( fastjet::PseudoJet( event[i].px(),
						event[i].py(), event[i].pz(), event[i].e() ) );
      }
    }
    
    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    vector <fastjet::PseudoJet> inclusiveAll, sortedAll;
    inclusiveAll = clustSeq.inclusive_jets(1.0);
    sortedAll    = sorted_by_pt(inclusiveAll);
    
    double ja_max = 0.20, jb_1 = 0.01, jb_2 = 0.03;

    if(sortedJets.size() > 0) {
      for(int j = 0; j < sortedJets.size(); j++){

	double scale = 1.0, res = 0.0;
	Vec4 pj = Vec4(sortedJets[j].px(),sortedJets[j].py(),
		       sortedJets[j].pz(),sortedJets[j].e());

	// Smearing
	double energy = pj.e(); 
	if(abs(sortedJets[j].rap()) < 2.0)
	  scale += jb_1;
	else
	  scale += jb_2;

	if(energy < 50.0)
	  res = ja_max * energy;
	else
	  res = (max(-0.002 * energy, -0.2) + 0.3) * energy;

	scale *= (1.0 + res * rndm.gauss()/energy);

	Vec4 pnew = scale * pj;

	bool placed = false;
	for(unsigned long j2 = 0; j2 < jets.size(); j2++){
	  if(jets.at(j2).pT() < pnew.pT()) {
	    jets.insert(jets.begin()+j2,pnew);
	    placed = true;
	    break;
	  }
	}
	if(!placed) jets.push_back(pnew);

	MET -= pnew;
      }
    }

    Vec4 unclustered(0.0, 0.0, 0.0, 0.0);
    for (int ij = sortedJets.size(); ij < sortedAll.size(); ij++){
      Vec4 pj = Vec4(sortedAll[ij].px(),sortedAll[ij].py(),
		     sortedAll[ij].pz(),sortedAll[ij].e());
      unclustered -= pj;
    }

    // TODO: smear the unclustered contribution
    
    MET += unclustered;
    
    
    // Naive Missing energy
    METnaive = Vec4(0.,0.,0.,0.);
    for (int j = 0; j < event.size(); ++j) {
      if (event[j].isFinal() && abs(event[j].eta()) < 4.5 && event[j].isVisible()){
	METnaive -= event[j].p();
      }
    }
    
    return true;
    
    // End ToyDetector::getObjects()
}

void ToyDetector::printObjects(){

    if(jets.size() > 0){
      cout << "Jets (pT, eta, phi)"<< endl;
      for(unsigned long i=0; i<jets.size(); i++){
	Vec4 obj = jets.at(i); 
	cout << obj.pT() <<"  "
	     << obj.eta() <<"  "
	     << obj.phi() <<"  "
	     << endl;
      }
    }

    if(electrons.size() > 0){
      cout << "Electrons (pT, eta, phi)"<< endl;
      for(unsigned long i=0; i<electrons.size(); i++){
    	Vec4 obj = electrons.at(i); 
    	cout << obj.pT() <<"  "
    	     << obj.eta() <<"  "
    	     << obj.phi() <<"  "
    	     << endl;
      }
    }
    if(muons.size() > 0){
      cout << "Muons (pT, eta, phi)"<< endl;
      for(unsigned long i=0; i<muons.size(); i++){
    	Vec4 obj = muons.at(i); 
    	cout << obj.pT() <<"  "
    	     << obj.eta() <<"  "
    	     << obj.phi() <<"  "
    	     << endl;
      }
    }

    if(MET.pT() > 0.0) {
      cout << "Missing energy (pT, phi)"<< endl;
      cout << MET.pT() <<"  "
	   << MET.phi() <<"  "
	   << endl;
      cout << "Naive Missing energy (pT, phi)"<< endl;
      cout << METnaive.pT() <<"  "
	   << METnaive.phi() <<"  "
	   << endl;
    }
    cout << "----------"<<endl;

    // End ToyDetector::printObjects()
}


// End namespace Pythia 8
}



// void ToyDetector::sortByPT(vector<Vec4> vect){

//     for(int i=0; i<vect.size()-1; i++){
//       for(int j=i+1; j<vect.size(); j++){
// 	Vec4 obj1 = vect.at(i); 
// 	Vec4 obj2 = vect.at(j); 
// 	if(obj1.pT() < obj2.pT()) {
	  
// 	}
//       }
      
//     }
// }
