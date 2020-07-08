## Event Acceptance x Efficiency
Event_AE = float(N_event_at_Event_selection)/float(LoadEvents)
## Number of charginos at E-S level
N_tracklets_at_ES = float(N_charginos_at_Event_selection)

## Loop over tau array
for tau in range(0,len(tau_array)):
   ## Tracklet Acceptance x Efficiency
   Tracklet_AE[tau] = N_traclets_rec[tau]/N_tracklets_at_ES
   ## Total Acceptance x Efficiency
   EA_EE_TA_TE[tau] = Event_AE*Tracklet_AE[tau]
   ## convert from pb to fb
   from_pb_to_fb = 1000
   ##cross section of tracklets with PT > 100 GeV in (fb)
   xs_Tracklet_PT100[tau] = xs_PT100_rec[tau]*from_pb_to_fb
   ##Number of expected tracklets events with PT > 100 GeV
   N_Exp_Tracklets_PT100[tau] = xs_Tracklet_PT100[tau]*Lum
   ## Number of MC tracklets at reconstruction level
   MCN_tracklets_rec[tau] = xs_Tracklet_PT100[tau]/(weight_xs*from_pb_to_fb)
   fAccEff.write('  {:.3f}           {:.3f}         {:.3f}      {:.4e}      {:.4e}      {:.4e}    {:.3e}    {:.2e}    {:.2e} \n '.format(Chargino_mass, Neutralino_mass, tau_array[tau], Event_AE, Tracklet_AE[tau], EA_EE_TA_TE[tau], xs_Tracklet_PT100[tau], N_Exp_Tracklets_PT100[tau], MCN_tracklets_rec[tau]))
