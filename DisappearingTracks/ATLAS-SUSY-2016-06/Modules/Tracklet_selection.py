## Tracklet PT efficiency (fix for all events)
TP = 0.57
## Total weight of the event
weight_event = weight_xs*(weight_11*TP + (weight_21*TP - weight_22*TP**2))
#probability for an event with N charginos to have at least one reconstructed tracklet
weight_HM = weight_11 + weight_21 - weight_22
###################################################
if(tau_array[t]):
  N_traclets_rec[t] += weight_HM
  if(MAX_PT>100):
	xs_PT100_rec[t] += weight_event
	MCevent[t] += 1.0
###################################################
