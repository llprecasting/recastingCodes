##################### Open Efficiency Map ########################
mapFilePath = Modules_dir+'DisappearingTrack2016-TrackAcceptanceEfficiency.root'
ElectroweakHistName = 'ElectroweakEfficiency'
###### Read efficiency map ############
mapFile = ROOT.TFile(mapFilePath)
tmap = mapFile.Get(ElectroweakHistName)
#######################################

##### Open variables file #######
exec(open(Modules_dir+'Internal_variables.py').read())
print'---- LLP Recast ---- '
print'Reading root file = ', rootFile

########### Read root file ############
froot = ROOT.TFile(rootFilePath,"READ")
tree = froot.Get("Delphes")

# Get the number of events of the root file
LoadEvents = tree.GetEntriesFast()
print"LoadEvents =", LoadEvents

##########################################################
##### main loop of the program, loop over the events #####
##########################################################
for jentry in xrange(LoadEvents):
    # Copy next entry into memory and verify.
    nb = tree.GetEntry( jentry )
    if nb<=0:
        continue
#### Counter of event, print a message every 1000 events ####
    if (m_eventCounter % 1000) == 0:
        print "Event number =", m_eventCounter
    m_eventCounter += 1

############### Event Selection-Level ########################
    weight_xs = tree.Event.At(0).Weight*kfac*1e9 # weight of each event
########## fill branches for Event selection #######
    OffMET	     = copy.deepcopy(tree.MissingET)
    jets50       = copy.deepcopy(tree.Jet)

    ########### Trigger on MissingET #################
    if( tree.MissingET.At(0).MET > 110 ):

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
        if( (electrons.GetEntries() == 0 or muons.GetEntries() == 0) ):

        ############## Event Selection ##################
            if( jets.GetEntries() != 0 ):
                if( jets.At(0) == None ): continue
                if( jets.At(0).PT >= 140 ): #leading jet
                    if( OffMET.At(0).MET >= 140 ):
                        if( dPhiMin(jets50,OffMET,4) >=1 ):
                            # N of events passing the Event Selection level
                            N_event_at_Event_selection += 1

                            ############## Tracklet Selection-Level ################
                            #Select charginos and neutralinos from Particle branch
                            charginos    = selectChargino(tree.Particle, 62, PID_chargino)
                            neutralinos  = selectNeutralino(tree.Particle, 1, PID_neutralino)
                            ## Filter: cuts on PT, Eta ##
                            charginos   = filterPhaseSpace(charginos, 5, 2.5)

                            # Chargino and neutralino mass
                            for c,n in zip(charginos, neutralinos):
                                Neutralino_mass = n.Mass
                                Chargino_mass = c.Mass

                            ## number of charginos per event
                            nC = charginos.GetEntries()
                            #N of charginos passing the E-S level
                            N_charginos_at_Event_selection += nC

                            ## loop over the lifetime
                            for t in range(0,len(tau_array)):
                                tau = tau_array[t]*1e-9 #convert to nano-seconds
                                #loop to increase statistic
                                for l in range(0, chargino_loop):
                                ## events with 1 chargino
                                 if(nC == 1):
                                    for c in charginos:
                                        R = CDF(tau, c.PT, c.Mass) ## radial decay length
                                        weight_HM = tmap.GetBinContent(tmap.FindBin(c.Eta,R))/chargino_loop ## weight Heatmap
                                        weight_11 = weight_HM # weight for events with 1 chargino
                                        weight_21 = 0		  # weight for events with 2 chargino selecting 1
                                        weight_22 = 0	      # weight for events with 2 chargino selecting 2
                                        MAX_PT = c.PT
                                        exec(open(Modules_dir+'Tracklet_selection.py').read())

                                ## events with 2 charginos
                                 if(nC == 2):
                                    cPT2 = []  # array to store the PT of each chargino
                                    R2 = []    # array to store the radial decay length of each chargino
                                    w_HM2 = [] # array to store the weight of Heatmap of each chargino
                                    for c in charginos:
                                        cPT2.append(c.PT)
                                        R = CDF(tau, c.PT, c.Mass) ## radial decay length
                                        ###############################################
                                        R2.append(R)
                                        w_HM2.append(tmap.GetBinContent(tmap.FindBin(c.Eta,R))/chargino_loop)
                                    ## Probability of selecting tracklets in case of 2 charginos
                                    weight_11 = 0                   # weight for events with 1 chargino
                                    weight_21 = w_HM2[0] + w_HM2[1] # weight for events with 2 chargino selecting 1
                                    weight_22 = w_HM2[0]*w_HM2[1]   # weight for events with 2 chargino selecting 2
                                    MAX_PT = max(cPT2[0],cPT2[1])   # selecting the chargino with higher PT
                                    exec(open(Modules_dir+'Tracklet_selection.py').read())

########### end of main loop ###########
