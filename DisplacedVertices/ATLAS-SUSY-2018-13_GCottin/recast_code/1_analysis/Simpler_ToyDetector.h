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

#include <Pythia8/Pythia.h>
#include <Pythia8/Event.h>
#include <vector>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

namespace Pythia8{

class ToyDetector{

 public:

  // Constructors
  ToyDetector(){ }

  bool getObjects(const Event& event);
  void printObjects();
  void clearAll();

  vector <Vec4> leptons, jets, electrons, muons;
  Vec4 MET, METnaive;

  // Fastjet analysis - select algorithm and parameters
  static double Rparam;
  static fastjet::Strategy strategy;
  static fastjet::RecombinationScheme recombScheme;
  static fastjet::JetDefinition *jetDef;

  static bool DEBUG;
  
 private:
  
  // Fastjet input
  vector <fastjet::PseudoJet> fjInputs;

  Rndm rndm;

  // End ToyDetector
};

// End namespace Pythia 8
}
