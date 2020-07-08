# Variables of the analysis:

## Definition of variables to use
Chargino_mass = 0		    ## Chargino mass
Neutralino_mass = 0		    ## Neutralino mass
m_eventCounter = 0  		## Counter of events
chargino_loop = 100		    ## Chargino loop to increase statistic
N_event_at_Event_selection = 0.     ##Number of events at Event-Selection level
N_charginos_at_Event_selection = 0. ##Number of charginos at Event-Selection level
N_traclets_rec = np.array([0.]*len(tau_array)) ## Number of traclets at reconstruction level
Tracklet_AE = np.array([0.]*len(tau_array))    ## Tracklet Aceptance x Efficiency
EA_EE_TA_TE = np.array([0.]*len(tau_array))    ## Total Acceptance x Efficiency
xs_PT100_rec = np.array([0.]*len(tau_array))   ## Total cross section of tracklets with PT>100 GeV in pb
MCN_tracklets_rec = np.array([0.]*len(tau_array)) ## Number of MC tracklets at reconstruction level
xs_Tracklet_PT100 = np.array([0.]*len(tau_array)) ## Total cross section of tracklets with PT>100 GeV in fb
N_Exp_Tracklets_PT100 = np.array([0.]*len(tau_array)) ## Number of expected tracklets events with PT > 100 GeV
MCevent = np.array([0.]*len(tau_array))
