#!/usr/bin/env python

from __future__ import print_function
import ROOT
import copy, os, sys
from Functions import filterPhaseSpace, overlapRemovalInterval, overlapRemoval, OffMissingET
from Functions import dPhiMin, selectChargino, selectNeutralino, CDF
import numpy as np


class EventSelector(object):


    def __init__(self,tau_array,PID_chargino=1000024,PID_neutralino=1000022,
                    chargino_loop=100,kfactor=1.0,lum=36.1):

        self.PID_chargino = PID_chargino
        self.PID_neutralino = PID_neutralino
        self.chargino_loop = chargino_loop #Controls how many times the decay length distribution is sampled (for a given proper lifetime)
        self.tau_array = tau_array
        self.kfac = kfactor
        self.lum = lum

        self.loadEfficiencies()

    def loadEfficiencies(self):

        #Efficiency map file:
        effFile = './Modules/DisappearingTrack2016-TrackAcceptanceEfficiency.root'
        if not os.path.isfile(effFile):
            print('Could not find ATLAS efficiencies file: %s ' %effFile)
            sys.exit()

        self.effFile = os.path.abspath(effFile)

    def loadDelphesLib(self,delphesPath):

        if not os.path.isdir(delphesPath):
            print('Can not find folder %s' %delphesPath)
            sys.exit()
        libPath = os.path.abspath(os.path.join(delphesPath,'libDelphes.so'))
        if not os.path.isfile(libPath):
            print('Can not find %s' %libPath)
            sys.exit()
        includePath = os.path.abspath(os.path.join(delphesPath,'external'))
        if not os.path.isdir(includePath):
            print('Can not find %s' %includePath)
            sys.exit()
        if not 'ROOT_INCLUDE_PATH' in os.environ or os.environ['ROOT_INCLUDE_PATH'] != includePath:
            print('Set the ROOT_INCLUDE_PATH enviroment variable to %s'  %(includePath))
            print('(i.e. export ROOT_INCLUDE_PATH=%s)' %(includePath))
            sys.exit()

        ROOT.gSystem.Load(libPath.replace('libDelphes.so','libDelphes'))

    def loadRootFile(self,rootFile):

        self.rootFile = os.path.abspath(rootFile)

    def resetVars(self):

        #Store the number of tracks for each lifetime:
        self.N_tracklets_rec = np.array([0.]*len(self.tau_array))
        #Store the event weights for each lifetime:
        self.xs_PT100_rec = np.array([0.]*len(self.tau_array))
        self.N_tracklets_rec = np.array([0.]*len(self.tau_array))
        #Track the number of events which passes the selection:
        self.MCevents = 0
        #Track the total cross-section of the events read:
        self.totalxsec = 0
        #Track the total cross-section of the events selected:
        self.selectxsec = 0
        #Store the current event weight
        self.weight_xs = 0
        #Track the number of events read from file
        self.EventsRead = 0
        #Track total number of charginos:
        self.N_charginos_at_Event_selection = 0
        #Track number of events after pre-selection:
        self.N_event_at_Event_selection = 0

        #Store neutralino mass
        self.neutralinoMass = None
        #Store chargino mass
        self.charginoMass = None

    def selectEvents(self,nmax=-1):

        self.resetVars()

        froot = ROOT.TFile(self.rootFile,"READ")
        tree = froot.Get("Delphes")

        ElectroweakHistName = 'ElectroweakEfficiency'
        ###### Read efficiency map ############
        mapFile = ROOT.TFile(os.path.abspath(self.effFile),"READ")
        effMap = mapFile.Get(ElectroweakHistName)

        #Total number of events:
        LoadEvents = tree.GetEntriesFast()

        #If nmax is positive just read up to nmax events
        if nmax <= 0:
            nEvt = LoadEvents
        else:
            nEvt = min(nmax,LoadEvents)

        print('Reading %i events' %nEvt)
        #Event loop:
        for ievent in range(nEvt):
            #### Counter of event, print a message every 1000 events ####
            if (ievent % 1000) == 0:
                print( "Reading event", ievent)
            # Copy next entry into memory and verify.
            nb = tree.GetEntry( ievent )
            if nb<=0:
                return
            self.EventsRead += 1
            preSel = self.preSelectEvent(tree)
            #If failed pre-selection, continue
            if not preSel:
                continue
            weights = self.computeEventWeights(tree,effMap)
            self.xs_PT100_rec += weights[0]
            self.N_tracklets_rec += weights[1]

    def preSelectEvent(self,tree):

        ########## fill branches for Event selection #######
        OffMET	     = copy.deepcopy(tree.MissingET)
        jets50       = copy.deepcopy(tree.Jet)

        ########### Trigger on MissingET #################
        if( tree.MissingET.At(0).MET < 110 ):
            return False

        ############## Event Reconstruction ###################
        ## prefilter: cuts on PT, Eta ##
        jets        = filterPhaseSpace(tree.Jet,       20, 2.8)
        jets50      = filterPhaseSpace(jets50,         50, 2.5)
        muons       = filterPhaseSpace(tree.Muon,      10, 2.7)
        electrons   = filterPhaseSpace(tree.Electron,  10,2.47)
        #charginos   = filterPhaseSpace(charginos,       5, 2.5)

        ## overlap removal: candidate, neighbours, dR_min, dR_max ##
        muons     = overlapRemovalInterval(muons,     jets,   0.2, 0.4)
        electrons = overlapRemovalInterval(electrons, jets,   0.2, 0.4)

        ## overlap removal:          candidate,    neighbours,   dR ##
        jets      = overlapRemoval(jets,         electrons,    0.2)
        jets      = overlapRemoval(jets,             muons,    0.2)

        #Offline MissingET ##
        OffMET = OffMissingET(jets, electrons, muons, OffMET) # Offline Missing ET

        ## lepto veto ##
        if( (electrons.GetEntries() != 0 or muons.GetEntries() != 0) ):
            return False

        ############## Event Selection ##################
        if( jets.GetEntries() == 0 ):
            return False
        if( jets.At(0) == None ):
            return False
        if( jets.At(0).PT < 140 ): #leading jet
            return False
        if( OffMET.At(0).MET < 140 ):
            return False
        if( dPhiMin(jets50,OffMET,4) < 1 ):
            return False

        charginos   = selectChargino(tree.Particle, 62, self.PID_chargino)
        charginos   = filterPhaseSpace(charginos, 5, 2.5)
        # number of charginos per event
        nC = charginos.GetEntries()
        if (nC > 2 or nC == 0):
            if (nC > 2):
                print('Can not deal with events with %i charginos' %nC)
            return False

        return True

    def computeEventWeights(self,tree,effMap):

        kfac = self.kfac

        ############### Event Selection-Level ########################
        weight_xs = tree.Event.At(0).Weight*kfac*1e9 # weight of each event
        self.totalxsec += weight_xs
        self.weight_xs = weight_xs

        ############## Tracklet Selection-Level ################
        #Select charginos and neutralinos from Particle branch
        charginos    = selectChargino(tree.Particle, 62, self.PID_chargino)
        neutralinos  = selectNeutralino(tree.Particle, 1, self.PID_neutralino)
        ## Filter: cuts on PT, Eta ##
        charginos   = filterPhaseSpace(charginos, 5, 2.5)
        ## number of charginos per event
        nC = charginos.GetEntries()

        #Count total number of charginos
        self.N_charginos_at_Event_selection += nC
        #Count number of events after pre-selection
        self.N_event_at_Event_selection += 1

        #Store neutralino and chargino masses (only used in the output)
        if self.neutralinoMass is None:
            for n in neutralinos:
                self.neutralinoMass = n.Mass
                break
        if self.charginoMass is None:
            for c in charginos:
                self.charginoMass = c.Mass
                break

        MAX_PT = max([c.PT for c in charginos])
        if MAX_PT > 100:
            self.selectxsec += weight_xs
            self.MCevents += 1

        #Compute weights for all lifetimes
        event_weights = [0.]*len(self.tau_array)
        ntracklets = [0.]*len(self.tau_array)
        for i,tauns in enumerate(self.tau_array):
            tau = tauns*1e-9 #convert to seconds
            weight_event,weight_Track = self.trackletSection(charginos,tau,effMap)
            ntracklets[i] = weight_Track
            if MAX_PT > 100:
                event_weights[i] = weight_xs*weight_event

        return (event_weights,ntracklets)

    def trackletSection(self,charginos,tauns,effMap):
        ###Compute probability of selecting at least one track

        #First build zero efficiencies dictionary
        #(assuming there are at most 2 charginos)
        Effvalues = dict([[i,np.array([0.]*self.chargino_loop)] for i in range(2)])
        ##sample the radial decay length
        # and compute efficiencies for each value:
        for iR in range(self.chargino_loop):
            for i,c in enumerate(charginos):
                R = CDF(tauns, c.PT, c.Mass) #decay length
                #Efficiency:
                Effvalues[i][iR] = effMap.GetBinContent(effMap.FindBin(c.Eta,R))

        #Now compute the event weight:
        #(actually list of weights for each decay length)
        weight_1 = Effvalues[0]+Effvalues[1] # weight for selecting at least 1 chargino
        weight_2 = Effvalues[0]*Effvalues[1]  # weight for selecting both charginos

        ## Tracklet PT efficiency (fix for all events)
        TP = 0.57
        ## Total weight of the event
        weight_event = weight_1*TP -weight_2*TP**2
        #probability for an event with N charginos to have at least one reconstructed tracklet
        weight_Track = weight_1 - weight_2

        #Now compute the weight averaged over all decay lengths:
        weight_event = np.average(weight_event)
        weight_Track = np.average(weight_Track)

        return (weight_event,weight_Track)

    def writeResults(self,outputFile):

        fAccEff = open(outputFile,'w')
        fAccEff.write(' Char_Mass(GeV)   Neutr_Mass(GeV)   tau(ns)      EAxEE          TAxTE         EAxEExTAxTE    xs100(fb)    Ntr        MCev  \n')

        from_pb_to_fb = 1000
        Chargino_mass = self.charginoMass
        Neutralino_mass = self.neutralinoMass

        ## Event Acceptance x Efficiency
        Event_AE = float(self.N_event_at_Event_selection)/float(self.EventsRead)
        ## Number of charginos at E-S level
        N_tracklets_at_ES = float(self.N_charginos_at_Event_selection)

        #Event average weight:
        # weight_xs = self.totalxsec/self.EventsRead
        #Weight of las event:
        weight_xs = self.weight_xs

        ## Loop over tau array
        for i,tau in enumerate(self.tau_array):
            ## Tracklet Acceptance x Efficiency
            Tracklet_AE = self.N_tracklets_rec[i]/N_tracklets_at_ES
            ## Total Acceptance x Efficiency
            EA_EE_TA_TE = Event_AE*Tracklet_AE
            ##cross section of tracklets with PT > 100 GeV in (fb)
            xs_Tracklet_PT100 = self.xs_PT100_rec[i]*from_pb_to_fb
            ##Number of expected tracklets events with PT > 100 GeV
            N_Exp_Tracklets_PT100 = xs_Tracklet_PT100*self.lum
            ## Number of MC tracklets at reconstruction level
            MCN_tracklets_rec = xs_Tracklet_PT100/(weight_xs*from_pb_to_fb)
            fAccEff.write('  {:.3f}           {:.3f}         {:.3f}      {:.4e}      {:.4e}      {:.4e}    {:.3e}    {:.2e}    {:.2e} \n '.format(Chargino_mass, Neutralino_mass, tau, Event_AE, Tracklet_AE, EA_EE_TA_TE, xs_Tracklet_PT100, N_Exp_Tracklets_PT100, MCN_tracklets_rec))


        fAccEff.close()
