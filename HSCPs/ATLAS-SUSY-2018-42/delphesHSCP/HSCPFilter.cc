/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class HSCPFilter
 *
 *  Filter particles with a given electric charge
 *  and minimal mass and pT values. Based on HSCPFilter.
 *
 *  \author A. Lessa - UFABC, Sao Paulo
 *
 */

#include "modules/HSCPFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

HSCPFilter::HSCPFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------


HSCPFilter::~HSCPFilter()
{
}

//------------------------------------------------------------------------------

void HSCPFilter::Init()
{

  ExRootConfParam pdgs,charges;
  Int_t i;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fMassMin = GetDouble("massMin", 20.0);

  charges = GetParam("Charge"); // Filter multiple charges

  pdgs = GetParam("PdgCode"); // Allow to specify PDGs

  fInvert = GetBool("Invert", false); // Allow to invert the filter

  fCharges.clear();
  for(i = 0; i < charges.GetSize(); ++i)
  {
    fCharges.push_back(charges[i].GetInt());
  }

  // If charges have not been defined, set to +-1
  if (fCharges.size() == 0){
    fCharges.push_back(1);
    fCharges.push_back(-1);
  }

  // Store PDF codes (if given)
  fPdgCodes.clear();
  for(i = 0; i < pdgs.GetSize(); ++i)
  {
    fPdgCodes.push_back(pdgs[i].GetInt());
  }

  // import input array iterator
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "hscpParticles"));

  fDaughtersArray = ExportArray(GetString("DaughtersArray", "daughters"));

}


//------------------------------------------------------------------------------

void HSCPFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void HSCPFilter::Process()
{
  Candidate *candidate;
  Candidate *newCandidate;
  Int_t pdgCode;
  Int_t i;
  Int_t charge;
  Bool_t pass;
  Double_t pt;
  Double_t mass;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pt = candidateMomentum.Pt();
    mass = candidate->Mass;
    charge = candidate->Charge;

    if(pt < fPTMin) continue;
    if(mass < fMassMin) continue;
    pass = kTRUE;

    // Select particles with PDG matching the input PDG list
    if (fPdgCodes.size() > 0){
      if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end()) pass = kFALSE;
    }
    if(fInvert) pass = !pass;

    // Select particles according to Charges list
    if(find(fCharges.begin(), fCharges.end(), charge) == fCharges.end()) pass = kFALSE;

    if (!pass) continue;

    // Now make sure that the particle is at the last step
    // of the event

    Candidate *daughter = 0;
    if(candidate->D1 < 0) continue;
    if(candidate->D2 < candidate->D1) continue;
    if(candidate->D1 >= fInputArray->GetEntriesFast() || candidate->D2 >= fInputArray->GetEntriesFast())
    {
      throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
    }

    for(i = candidate->D1; i <= candidate->D2; ++i)
    {
      daughter = static_cast<Candidate *>(fInputArray->At(i));
      // If any of the daughter PID matches the candidate
      // PDG, the candidate is an intermediate step
      if (daughter->PID == pdgCode){
        pass = kFALSE;
        break;
      };
    }

    if (!pass) continue;

    newCandidate = static_cast<Candidate *>(candidate->Clone());
    // Set daughter index according to the DaughtersArray:
    newCandidate->D1 = fDaughtersArray->GetEntries();
    // Store the daughters
    for(i = candidate->D1; i <= candidate->D2; ++i)
    {
      daughter = static_cast<Candidate *>(fInputArray->At(i)->Clone());
      // Redefined the daughters mother acconding to output array
      daughter->M1 = fOutputArray->GetEntries();
      daughter->M2 = fOutputArray->GetEntries();
      fDaughtersArray->Add(daughter);
    }
    newCandidate->D2 = fDaughtersArray->GetEntries()-1;
    fOutputArray->Add(newCandidate);
  }
}


