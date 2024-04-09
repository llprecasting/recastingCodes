/*
Recaasting code for:

Search for long-lived, massive particles in events with displaced vertices and multiple jets in pp collisions at \sqrt{s} = 13 TeV with the ATLAS detector
arXiv: 2301.13866

Kingman Cheung <cheung@phys.nthu.edu.tw>, Fei-Tung Chung <feitung.chung@gapp.nthu.edu.tw>, Giovanna Cottin <gfcottin@uc.cl>, Zeren Simon Wang <wzs@mx.nthu.edu.tw>

*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>

#include "Pythia8/Pythia.h"
#include "Simpler_ToyDetector.h"
#include "efficiency.h"

using namespace Pythia8;

int main(int argc, char* argv[])
{

	if (argc != 2) {
        	std::cout << "./recast_2301_13866 pathToCMND" << std::endl;
        	std::cout << "   - pathToCMND is the path to the Pythia input cmnd card" << std::endl;
        	exit(1);
    	}
    

	Pythia pythia;
	pythia.settings.addMode("LHEFInputs:nSubruns", 0, true, false, 0, 100); //allows LHEFInputs:nSubruns
	Event& event = pythia.event;

	std::string pathToCMND = argv[1];
	pythia.readFile(pathToCMND);

	int nEvent = pythia.mode("Main:numberOfEvents");
	int nAbort = pythia.mode("Main:timesAllowErrors");
	int iAbort = 0;

	double nHighPTJetSelection = 0.;
	double nFiducialHighPT = 0.;
	double nTransverseDistanceHighPT = 0.;
	double nTransverseImpactParameterHighPT = 0.;
	double nSelectedDecayProductsHighPT = 0.;
	double nInvariantMassHighPT = 0.;

	double nTracklessJetSelection = 0.;
	double nFiducialTrackless = 0.;
	double nTransverseDistanceTrackless = 0.;
	double nTransverseImpactParameterTrackless = 0.;
	double nSelectedDecayProductsTrackless = 0.;
	double nInvariantMassTrackless = 0.;

	double nHighPTAcc = 0.;
	double nTracklessAcc = 0.;
	double nHighPTAccEff = 0.;
	double nTracklessAccEff = 0.;

	ToyDetector detector;


	//Extract number of events and max number of jets in merging.

	int nMerge = pythia.mode("LHEFInputs:nSubruns") - 1; //pythia.mode("Merging:nJetMax");
	//cout << "nMerge: " << nMerge << endl;

	//Merged total cross section, summed over subruns.
	double sigmaTotal = 0.;
	double allevents = 0;

	//Loop over subruns with varying number of jets.
	for (int iMerge = 0; iMerge <= nMerge; ++iMerge)
	{
		double sigmaSample = 0.;

		//Read in name of LHE file for current subrun and initialize.
		pythia.readFile(pathToCMND, iMerge);

		pythia.init();


		//Event loop begin
		for (int iEvent = 0; iEvent < nEvent; iEvent++)
		{
			if (!pythia.next())
			{
				if (pythia.info.atEndOfFile())
				{
					std::cout << "Reached end of Les Houches Event File!" << std::endl;
					break;
				}
				if (++iAbort < nAbort)
					continue;
				std::cout << "Too many errors!" << std::endl;
				return 1;
			}


			//Get CKKWL weight of current event.
			double evtweight = pythia.info.weight();
			//cout  << scientific << setprecision(4)<< "evtweight: " << evtweight << endl;
			double weight = pythia.info.mergingWeight();

			//cout  << scientific << setprecision(4)<< "first weight: " << weight << endl;

			weight *= evtweight * 1e-9; //pb convereted to mb
			//cout  << scientific << setprecision(4)<< "second weight: " << weight << endl;
			sigmaSample += weight;
			allevents += weight;

			//Do nothing for vanishing weight (event record might not be filled)
			if (weight == 0)
				continue;


			//Jet selection

			detector.getObjects(event);
			double sumPT = 0.;
			for (int iJet = 0; iJet < detector.jets.size(); iJet++)
			{
				sumPT += detector.jets[iJet].pT();
			}

			//Find mothers

			std::vector<int> motherIndex; //Tracks of mothers
			motherIndex.clear();
			for (int iTrack = 0; iTrack < event.size(); iTrack++)
			{
				if ((event[iTrack].idAbs() == 1000022) || (event[iTrack].idAbs() == 1000023) || (event[iTrack].idAbs() == 1000024))
				//if (event[iTrack].idAbs() == 9900035)
					motherIndex.push_back(iTrack);
			}

			//Jet selection

			std::vector<int> matchedJetIndex;
			std::vector<int> matchedTrackIndex;
			std::vector<double> displacedJetsPT;
			matchedTrackIndex.clear();
			matchedJetIndex.clear();
			displacedJetsPT.clear();
			for (int iMother : motherIndex)
			{
				std::vector<int> daughterIndex = event[iMother].daughterListRecursive();
				std::sort(daughterIndex.begin(), daughterIndex.end());
				daughterIndex.erase(std::unique(daughterIndex.begin(), daughterIndex.end()), daughterIndex.end());
				for (int jJet = 0; jJet < detector.jets.size(); jJet++)
				{
					double minDeltaR = std::numeric_limits<double>::max();
					int matchedTrack = -1;
					for (int kTrack : daughterIndex)
					{
						if (!event[kTrack].isCharged() || !event[kTrack].isFinal())
							continue;
						double deltaEta = detector.jets[jJet].eta() - event[kTrack].eta();
						double deltaPhi = detector.jets[jJet].phi() - event[kTrack].phi();
						if (std::abs(deltaPhi) > M_PI)
							deltaPhi = 2. * M_PI - std::abs(deltaPhi);
						double deltaR = std::sqrt(std::pow(deltaEta, 2.) + std::pow(deltaPhi, 2.));
						if ((deltaR < 0.3) && (deltaR < minDeltaR))
						{
							minDeltaR = deltaR;
							matchedTrack = kTrack;
						}
					}
					if (matchedTrack >= 0)
					{
						matchedJetIndex.push_back(jJet);
						matchedTrackIndex.push_back(matchedTrack);
						//std::cout << event[iMother].xDec() << event[iMother].yDec() << event[iMother].zDec() << " " << event[matchedTrack].xProd() << event[matchedTrack].yProd() << event[matchedTrack].zProd() << std::endl;
					}
				}
			}

			for (int i = 0; i < matchedJetIndex.size(); i++)
			{
				int jetIndex = matchedJetIndex[i];
				int trackIndex = matchedTrackIndex[i];
				if (std::abs(detector.jets[jetIndex].eta()) > 2.5)
					continue;
				double rVertex = std::sqrt(std::pow(event[event[trackIndex].mother1()].xDec(), 2.) + std::pow(event[event[trackIndex].mother1()].yDec(), 2.));
				//std::cout << event[event[trackIndex].mother1()].xDec() << event[event[trackIndex].mother1()].yDec() << event[event[trackIndex].mother1()].zDec() << " "
				//<< event[trackIndex].xProd() << event[trackIndex].yProd() << event[trackIndex].zProd() << std::endl;
				if (rVertex > 3870.)
					continue;
				displacedJetsPT.push_back(detector.jets[jetIndex].pT());
			}

			bool passHighPTJetSelection = false;
			bool passTracklessJetSelection = false;

			if ((detector.jets.size() >= 4) && (detector.jets[3].pT() > 250.))
				passHighPTJetSelection = true;
			if ((detector.jets.size() >= 5) && (detector.jets[4].pT() > 195.))
				passHighPTJetSelection = true;
			if ((detector.jets.size() >= 6) && (detector.jets[5].pT() > 116.))
				passHighPTJetSelection = true;
			if ((detector.jets.size() >= 7) && (detector.jets[6].pT() > 90.))
				passHighPTJetSelection = true;

			if ((detector.jets.size() >= 4) && (detector.jets[3].pT() > 137.))
				passTracklessJetSelection = true;
			if ((detector.jets.size() >= 5) && (detector.jets[4].pT() > 101.))
				passTracklessJetSelection = true;
			if ((detector.jets.size() >= 6) && (detector.jets[5].pT() > 83.))
				passTracklessJetSelection = true;
			if ((detector.jets.size() >= 7) && (detector.jets[6].pT() > 55.))
				passTracklessJetSelection = true;

			bool passTracklessAdditional = false;

			std::sort(displacedJetsPT.begin(), displacedJetsPT.end(), std::greater<double>());
			if ((displacedJetsPT.size() >= 1) && (displacedJetsPT[0] > 70.))
				passTracklessAdditional = true;
			if ((displacedJetsPT.size() >= 2) && (displacedJetsPT[1] > 50.))
				passTracklessAdditional = true;


			bool passJetSelection = passHighPTJetSelection || passTracklessJetSelection && passTracklessAdditional;
			if (!passJetSelection)
				continue; //run faster!

			//Look through vertices and apply selection

			bool passFiducial = false;
			bool passTransverseDistance = false;
			bool passTransverseImpactParameter = false;
			bool passSelectedDecayProducts = false;
			bool passInvariantMass = false;
			bool passVertexAcc = false;
			double rMax = 0.; //Farthest truth LLP Decay
			double vEff = 1.; //Vertex efficiency of the event
			//Vertex loop begin
			for (int iMother : motherIndex)
			{
				double xVertex = event[iMother].xDec();
				double yVertex = event[iMother].yDec();
				double zVertex = event[iMother].zDec();
				double rVertex = std::sqrt(std::pow(xVertex, 2.) + std::pow(yVertex, 2.));
				rMax = std::max(rMax, rVertex);

				std::vector<int> daughterIndex = event[iMother].daughterListRecursive();
				std::sort(daughterIndex.begin(), daughterIndex.end());
				daughterIndex.erase(std::unique(daughterIndex.begin(), daughterIndex.end()), daughterIndex.end());
				int trackCount = 0;
				int impactedCount = 0;
				int selectedCount = 0;
				Vec4 totalMomentumAssumingPion(0., 0., 0., 0.);
				double const mChargedPion = 0.13957039;
				for (int jTrack : daughterIndex)
				{
					if (!event[jTrack].isCharged() || !event[jTrack].isFinal())
						continue;

					double rTrack = std::sqrt(std::pow(event[jTrack].xProd(), 2.) + std::pow(event[jTrack].yProd(), 2.));
					double zTrack = event[jTrack].zProd();
					double deltaPhi = event[jTrack].vProd().phi() - event[jTrack].phi();
					if (std::abs(deltaPhi) > M_PI)
						deltaPhi = 2. * M_PI - std::abs(deltaPhi);
					double transverseImpactParameter = rTrack * std::sin(deltaPhi);
					if (std::abs(transverseImpactParameter) > 2.)
						impactedCount++;

					double beta = event[jTrack].pT() / event[jTrack].e();
					double gamma = event[jTrack].e() / event[jTrack].m();
					double transverseDistance = beta * gamma * event[jTrack].tau();
					if ((transverseDistance > 520.) && (event[jTrack].pT() / std::abs(event[jTrack].charge()) > 1.))
						//if ((std::abs(event[jTrack].status()) == 91) && (event[jTrack].pT() / std::abs(event[jTrack].charge()) > 1.))
						selectedCount++;
					else
						continue; //Only selected decay products are used in calculating invariant mass of truth vertex

					double eAssumingPion = std::sqrt(event[jTrack].pAbs2() + std::pow(mChargedPion, 2.));
					totalMomentumAssumingPion += Vec4(event[jTrack].px(), event[jTrack].py(), event[jTrack].pz(), eAssumingPion);
				}

				//Fiducial
				if ((rVertex < 300.) && (std::abs(zVertex) < 300.))
					passFiducial = true;
				else
					continue;

				//Transverse distance
				if (rVertex > 4.)
					passTransverseDistance = true;
				else
					continue;

				//Transverse impact parameter
				if (impactedCount >= 1)
					passTransverseImpactParameter = true;
				else
					continue;

				//Selected decay products
				if (selectedCount >= 5)
					passSelectedDecayProducts = true;
				else
					continue;

				//Invariant mass
				if (totalMomentumAssumingPion.mCalc() > 10.)
					passInvariantMass = true;
				else
					continue;

				passVertexAcc = true;
				vEff *= 1. - vertexEff(rVertex, totalMomentumAssumingPion.mCalc(), selectedCount);
			}
			//Vertex loop end
			vEff = 1. - vEff;

			if (passHighPTJetSelection)
			{
				nHighPTJetSelection += weight;

				if (passFiducial)
					nFiducialHighPT += weight;
				if (passTransverseDistance)
					nTransverseDistanceHighPT += weight;
				if (passTransverseImpactParameter)
					nTransverseImpactParameterHighPT += weight;
				if (passSelectedDecayProducts)
					nSelectedDecayProductsHighPT += weight;
				if (passInvariantMass)
					nInvariantMassHighPT += weight;

				if (passVertexAcc)
				{
					nHighPTAcc += weight;
					nHighPTAccEff += weight * eventEff_HighPT(rMax, sumPT) * vEff;
				}
			}
			if (passTracklessJetSelection && passTracklessAdditional)
			{
				nTracklessJetSelection += weight;

				if (passFiducial)
					nFiducialTrackless += weight;
				if (passTransverseDistance)
					nTransverseDistanceTrackless += weight;
				if (passTransverseImpactParameter)
					nTransverseImpactParameterTrackless += weight;
				if (passSelectedDecayProducts)
					nSelectedDecayProductsTrackless += weight;
				if (passInvariantMass)
					nInvariantMassTrackless += weight;

				if (passVertexAcc)
				{
					nTracklessAcc += weight;
					nTracklessAccEff += weight * eventEff_Trackless(rMax, sumPT) * vEff;
				}
			}
		}
		//Event loop end

		pythia.stat();
		
		//Sum up merged cross section of current run.
		//cout << "sigmaSample: " << sigmaSample << endl;
		sigmaSample *= 1 / double(pythia.info.nAccepted()); //pythia.info.sigmaGen() / double(pythia.info.nAccepted());
		//cout << "pythia.info.sigmaGen(): " << pythia.info.sigmaGen() << endl;
		//cout << "pythia.info.nAccepted(): " << pythia.info.nAccepted() << endl;
		sigmaTotal += sigmaSample;
		//cout << "sigmaTotal at iMerge = " << iMerge << " is " << sigmaTotal << " mb." << endl;

	}

	using std::cout;
	using std::endl;
	cout << std::left;
	cout << std::setprecision(6);
	cout << endl;
	cout << std::setw(50) << "All Events: " << allevents / sigmaTotal << endl;
	cout << endl;
	cout << std::setw(50) << "High-pT Jet Selection [%]: " << nHighPTJetSelection / allevents * 100 << endl;
	cout << std::setw(50) << "Event has >= 1 DV passing [%]: " << endl;
	cout << std::setw(50) << "High-pT Fiducial [%]: " << nFiducialHighPT / allevents * 100 << endl;
	cout << std::setw(50) << "High-pT R_vertex > 4 mm [%]: " << nTransverseDistanceHighPT / allevents * 100 << endl;
	cout << std::setw(50) << "High-pT Have track with |d_0| > 2 mm [%]: " << nTransverseImpactParameterHighPT / allevents * 100 << endl;
	cout << std::setw(50) << "High-pT Selected Decay Product >= 5 [%]: " << nSelectedDecayProductsHighPT / allevents * 100 << endl;
	cout << std::setw(50) << "High-pT Invariant Mass > 10 GeV [%]: " << nInvariantMassHighPT / allevents * 100 << endl;
	cout << endl;
	cout << std::setw(50) << "Trackless Jet Selection [%]: " << nTracklessJetSelection / allevents * 100 << endl;
	cout << std::setw(50) << "Event has >= 1 DV passing [%]: " << endl;
	cout << std::setw(50) << "Trackless Fiducial [%]: " << nFiducialTrackless / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless R_vertex > 4 mm [%]: " << nTransverseDistanceTrackless / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless Have track with |d_0| > 2 mm [%]: " << nTransverseImpactParameterTrackless / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless Selected Decay Product >= 5 [%]: " << nSelectedDecayProductsTrackless / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless Invariant Mass > 10 GeV [%]: " << nInvariantMassTrackless / allevents * 100 << endl;
	cout << endl;
	cout << std::setw(50) << "High-pT Acc: " << nHighPTAcc / sigmaTotal << endl;
	cout << std::setw(50) << "High-pT Acc x Eff: " << nHighPTAccEff / sigmaTotal << endl;
	cout << std::setw(50) << "Trackless Acc: " << nTracklessAcc / sigmaTotal << endl;
	cout << std::setw(50) << "Trackless Acc x Eff: " << nTracklessAccEff / sigmaTotal << endl;
	cout << endl;
	cout << std::setw(50) << "High-pT Acc [%]: " << nHighPTAcc / allevents * 100 << endl;
	cout << std::setw(50) << "High-pT Acc x Eff [%]: " << nHighPTAccEff / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless Acc [%]: " << nTracklessAcc / allevents * 100 << endl;
	cout << std::setw(50) << "Trackless Acc x Eff [%]: " << nTracklessAccEff / allevents * 100 << endl;
	cout << endl;
	cout << std::setw(50) << "The inclusive cross section after merging is: " << sigmaTotal << " mb" << endl;

	return 0;
}
