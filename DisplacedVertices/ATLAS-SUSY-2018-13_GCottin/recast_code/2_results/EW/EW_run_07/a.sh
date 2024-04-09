#!/bin/bash

processNames=("c1pc1m" "n1c1m" "n1c1p" "n1n2" "n2c1m" "n2c1p")
signalRegion="High-pT"

nTotal=0.0
nJetSelection=0.0
nFiducial=0.0
nTransverseDistance=0.0
nTransverseImpactParameter=0.0
nSelectedDecayProducts=0.0
nInvariantMass=0.0
nAcc=0.0
nAccEff=0.0

for processName in "${processNames[@]}"
do
	xs=$(grep "The inclusive cross section after merging is: " ./$processName.txt | awk '{print $(NF-1)}')
	allEvents=$(grep "All Events: " ./$processName.txt | awk '{print $NF}')
	nTotal=$(awk "BEGIN { print($nTotal + $xs * $allEvents) }")
	nJetSelection=$(awk "BEGIN { print($nJetSelection + $xs * $allEvents * $(grep "$signalRegion Jet Selection \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nFiducial=$(awk "BEGIN { print($nFiducial + $xs * $allEvents * $(grep "$signalRegion Fiducial \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nTransverseDistance=$(awk "BEGIN { print($nTransverseDistance + $xs * $allEvents * $(grep "$signalRegion R_vertex > 4 mm \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nTransverseImpactParameter=$(awk "BEGIN { print($nTransverseImpactParameter + $xs * $allEvents * $(grep "$signalRegion Have track with |d_0| > 2 mm \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	#nTransverseImpactParameter=$(awk "BEGIN { print($nTransverseImpactParameter + $xs * $allEvents * $(grep "$signalRegion Have track with \|d_0\| > 2 mm \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nSelectedDecayProducts=$(awk "BEGIN { print($nSelectedDecayProducts + $xs * $allEvents * $(grep "$signalRegion Selected Decay Product >= 5 \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nInvariantMass=$(awk "BEGIN { print($nInvariantMass + $xs * $allEvents * $(grep "$signalRegion Invariant Mass > 10 GeV \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nAcc=$(awk "BEGIN { print($nAcc + $xs * $allEvents * $(grep "$signalRegion Acc \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
	nAccEff=$(awk "BEGIN { print($nAccEff + $xs * $allEvents * $(grep "$signalRegion Acc x Eff \[\%\]: " ./$processName.txt | awk '{print $NF}')) }")
done

> ./output.txt
echo $(awk "BEGIN { print($nJetSelection / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nFiducial / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nTransverseDistance / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nTransverseImpactParameter / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nSelectedDecayProducts / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nInvariantMass / $nTotal) }") >> ./output.txt
echo >> ./output.txt
echo $(awk "BEGIN { print($nAcc / $nTotal) }") >> ./output.txt
echo $(awk "BEGIN { print($nAccEff / $nTotal) }") >> ./output.txt
