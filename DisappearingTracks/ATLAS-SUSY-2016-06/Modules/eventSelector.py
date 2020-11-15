#!/usr/bin/env python

from __future__ import print_function
import ROOT
import copy, os, sys
from Functions import filterPhaseSpace, overlapRemovalInterval, overlapRemoval, OffMissingET
from Functions import dPhiMin, selectChargino, selectNeutralino, CDF
import numpy as np


class EventSelector(object):


    def __init__(self, parser):


        self.outputFile = parser.get('options','outputFile')
        self.rootFile = parser.get('options','rootFile')
        self.lum = parser.getfloat('options','Luminosity')
        self.PID_chargino = parser.getint('options','PIDchargino')
        self.PID_neutralino = parser.getint('options','PIDneutralino')
        self.kfac = parser.getfloat('options','kfactor')
        self.chargino_loop = parser.getint('options','nloop')
        self.weight = parser.get('options','Weight')
        self.effFile = parser.get('options','effFile')
        tau_values = parser.get('options','tau_values')
        try:
            tau_array = eval(tau_values)
        except:
            print("Define a valid expression for the tau_values parameter. It should be a list of values (e.g. [0.01,0.03,...]) or a single value (e.g. 0.01)")
            sys.exit()
        if isinstance(tau_array,(int,float)):
            self.tau_array = np.array([float(tau_array)])
        elif isinstance(tau_array,list):
            self.tau_array = np.array(tau_array)
        else:
            print('Define a valid value for the tau_values parameter. It should be a list or a float')
            sys.exit()

        if self.weight == 'same':
           if not parser.has_option('options','InitXs'):
               print('If you choose Weight: same, you must provide InitXs')
               sys.exit()
           else:
               self.init_xs = parser.getfloat('options','InitXs')
        elif self.weight != 'root':
            print('Weight option is not properly defined. Options are: same or root')
            sys.exit()


        self.loadEfficiencies()

    def loadEfficiencies(self):

        #Efficiency map file:
        effFile = self.effFile
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

    def loadRootFile(self):

        self.rootFile = os.path.abspath(self.rootFile)
        if not os.path.isfile(self.rootFile):
                print('Could not find root file')
                sys.exit()

    def resetVars(self):

        #Store the event weights for each lifetime:
        self.xs_PT100_rec = np.array([0.]*len(self.tau_array))
        #Store the numver of tracklets after reconstructuin level for each lifetime:
        self.N_tracklets_rec = np.array([0.]*len(self.tau_array))
        #Store the number of MC events which passes the tracklet selection:
        self.MCevents = np.array([0.]*len(self.tau_array))
        #Store the initial cross-section of the sample:
        self.init_xsec = 0
        #Store the cross-section of the events selected:
        self.selectxsec = 0
        #Store the current event weight
        self.weight_xs = 0
        #Store the total number of events of the sample
        self.EventsRead = 0
        #Store total number of charginos of the sample:
        self.N_charginos_at_Event_selection = 0
        #Store the number of events after pre-selection:
        self.N_event_at_Event_selection = 0
        #Store neutralino mass
        self.neutralinoMass = None
        #Store chargino mass
        self.charginoMass = None

    def selectEvents(self,nmax=-1):

        self.resetVars()

        kfac = self.kfac

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

            if(self.weight=='same'):
                weight_xs = (self.init_xs*kfac)/nEvt # all events have the same weight
            if(self.weight=='root'):
                weight_xs = tree.Event.At(0).Weight*kfac*1e9 # events with different weights

            self.init_xsec += weight_xs

            preSel = self.preSelectEvent(tree)
            #If failed pre-selection, continue
            if not preSel:
                continue
            #Count number of events after pre-selection
       	    self.N_event_at_Event_selection += 1
            self.selectxsec += weight_xs

            weights = self.computeEventWeights(tree, effMap, weight_xs)

            self.xs_PT100_rec += weights[0]
            self.N_tracklets_rec += weights[1]
            self.MCevents += weights[2]

    def preSelectEvent(self,tree):

        ########## fill branches for Event selection #######
        OffMET	     = copy.deepcopy(tree.MissingET)
        jets50       = copy.deepcopy(tree.Jet)

        ########### Trigger on MissingET #################
        if( tree.MissingET.At(0).MET <= 110 ):
            return False

        ############## Event Reconstruction ###################
        ## prefilter: cuts on PT, Eta ##
        jets        = filterPhaseSpace(tree.Jet,       20, 2.8)
        jets50      = filterPhaseSpace(jets50,         50, 2.5)
        muons       = filterPhaseSpace(tree.Muon,      10, 2.7)
        electrons   = filterPhaseSpace(tree.Electron,  10,2.47)

        ## overlap removal: candidate, neighbours, dR_min, dR_max ##
        muons     = overlapRemovalInterval(muons,     jets,   0.2, 0.4)
        electrons = overlapRemovalInterval(electrons, jets,   0.2, 0.4)

        ## overlap removal:          candidate,    neighbours,   dR ##
        jets      = overlapRemoval(jets,         electrons,    0.2)
        jets      = overlapRemoval(jets,             muons,    0.2)

        #Offline MissingET ##
        OffMET = OffMissingET(jets, electrons, muons, OffMET)

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

        return True

    def computeEventWeights(self,tree,effMap,weight_xs):

        ############## Tracklet Selection-Level ################
        #Select charginos and neutralinos from Particle branch
        charginos    = selectChargino(tree.Particle, 62, self.PID_chargino)
        neutralinos  = selectNeutralino(tree.Particle, 1, self.PID_neutralino)
        ## Filter: cuts on PT, Eta ##
        charginos   = filterPhaseSpace(charginos, 5, 2.5)
        ## number of charginos per event
        nC = charginos.GetEntries()

        if(nC == 0):
            return (0,0,0)

        #Count total number of charginos
        self.N_charginos_at_Event_selection += nC

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

        #Compute weights, num tracklets, num MC events and TAxTExTP for all lifetimes
        event_weights = [0.]*len(self.tau_array)
        ntracklets = [0.]*len(self.tau_array)
        nMCevents = [0.]*len(self.tau_array)
        for i,tauns in enumerate(self.tau_array): # loop sobre taus
            tau = tauns*1e-9 #convert to seconds
            weight_event, weight_Track, MC_event = self.trackletSection(charginos,tau,effMap)
            ntracklets[i] = weight_Track

            if MAX_PT > 100:
                event_weights[i] = weight_xs*weight_event
                if(weight_event != 0):
                    nMCevents[i] += MC_event

        return (event_weights,ntracklets,nMCevents)

    def trackletSection(self,charginos,tauns,effMap):
        ###Compute probability of selecting at least one track

        #First build zero efficiencies dictionary
        #(assuming there are at most 2 charginos)
        Effvalues = dict([[i,np.array([0.]*self.chargino_loop)] for i in range(2)])
        MCvalues = dict([[i,np.array([0.]*self.chargino_loop)] for i in range(2)])
        ##sample the radial decay length and compute efficiencies for each value:
        for iR in range(self.chargino_loop):
            for i,c in enumerate(charginos):
                R = CDF(tauns, c.PT, c.Mass) #decay length
                #Efficiency:
                HMeff = effMap.GetBinContent(effMap.FindBin(c.Eta,R))
                Effvalues[i][iR] = HMeff
                if(HMeff != 0.0):
                    MCvalues[i][iR] += 1

        #Now compute the event weight:
        #(actually list of weights for each decay length)
        weight_1 = Effvalues[0]+Effvalues[1] # weight for selecting at least 1 chargino
        weight_2 = Effvalues[0]*Effvalues[1]  # weight for selecting both charginos

        #(actually list of MC events for each decay length)
        MCevent = MCvalues[0]+MCvalues[1]

        ## Tracklet PT efficiency (fix for all events)
        TP = 0.57
        ## Total weight of the event
        weight_event = weight_1*TP - weight_2*TP**2
        # Probability for an event with N charginos to have at least one reconstructed tracklet
        weight_Track = weight_1 - weight_2
        # Now compute the weight averaged over all decay lengths:
        weight_event = np.average(weight_event)
        weight_Track = np.average(weight_Track)
        # Compute the total number of MC for each decay length
        MCevent_sum = np.sum(MCevent)

        return (weight_event, weight_Track, MCevent_sum)

    def writeResults(self):

        fAccEff = open(self.outputFile,'w')
        fAccEff.write(' C_Mass(GeV)   N_Mass(GeV)    tau(ns)    Init_xs(fb)     EAxEE          TAxTE         EAxEExTAxTE      total_eff     xs100(fb)     xslim(fb)         Ntr          r       MCev \n')

        from_pb_to_fb = 1000
        Chargino_mass = self.charginoMass
        Neutralino_mass = self.neutralinoMass


        ## Event Acceptance x Efficiency
        Event_AE = float(self.N_event_at_Event_selection)/float(self.EventsRead)
        ## Number of charginos at E-S level
        N_tracklets_at_ES = float(self.N_charginos_at_Event_selection)
        #Initial cross section of the sample:
        init_xsec_fb = self.init_xsec*from_pb_to_fb
        #Upper limit on the model-independent visible cross-section at 95% CL in fb. Table 4, article 1712.02118
        xs95lim = 0.22

        if self.N_event_at_Event_selection == 0:
            print('No events passed the selection.\n')
            fAccEff.write('No events passed the selection.\n')
            fAccEff.close()
            return
        elif N_tracklets_at_ES == 0:
            print('No charginos could be reconstructed.\n')
            fAccEff.write('No charginos could be reconstructed.\n')
            fAccEff.close()
            return


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
            N_MC_events_rec = self.MCevents[i]
            ## Total efficiency
            total_eff = xs_Tracklet_PT100/init_xsec_fb

            if(total_eff == 0.0):
                print('Warning: Efficiency of point Mass={:.2f}, tau={:.3f} is zero.\n More events are needed. Increase nloop variable.'.format(Chargino_mass,self.tau_array[i]))
                ## cross section limit
                xslim = np.nan
            else:
                xslim = xs95lim/(total_eff)

            ## r value
            r = N_Exp_Tracklets_PT100/(xs95lim*self.lum)
            ## cross section of events passing the Event Selection
            xs_select = self.selectxsec*from_pb_to_fb


            fAccEff.write('  {: .3f}     {:.3f}         {:.3f}      {:.3f}      {:.4e}      {:.4e}      {:.4e}      {:.4e}    {:.3e}     {:.3e}       {:.2e}    {:.2e}    {:.0f} \n '.format(Chargino_mass, Neutralino_mass, tau, init_xsec_fb, Event_AE, Tracklet_AE, EA_EE_TA_TE, total_eff, xs_Tracklet_PT100, xslim, N_Exp_Tracklets_PT100, r, N_MC_events_rec))

        fAccEff.close()
