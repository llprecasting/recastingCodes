# This code write script that MG and Pythia will read to generate events
import os
import sys

mode = sys.argv[1]
mass_Alp = float(sys.argv[2])
nevent = int(sys.argv[3])
slug = f'mALP{mass_Alp:.6g}_{mode}' 
InDir = sys.argv[4]
OutDir = sys.argv[5]


f = open(f"{OutDir}/script_{slug}.txt", 'w') # Creation of the Script to follow for MG.
#f.write(f"import model {InDir}/HAHM_MG5model_v3/HAHM_variableMW_v3_UFO \n") # Import the model.
#f.write(f"define f = u c d s u~ c~ d~ s~ b b~ e+ e- mu+ mu- ta+ ta- t t~ \n") # Define a fermion.
#f.write(f"generate g g > h HIG=1 HIW=0 QED=0 QCD=0, (h > h2 h2, h2 > f f) \n") # Generate the process.
f.write(f"import model sm \n")
f.write(f"import model {InDir}/ALP_linear_UFO_WIDTH \n")
f.write(f"define p = g u c d s u~ c~ d~ s~ \n")
f.write(f"define j = u c d s b u~ c~ d~ s~ b~ \n")
f.write(f"define l+ = e+ mu+ \n")
f.write(f"define l- = e- mu- \n")
f.write(f"define vl = ve vm \n")
f.write(f"define vl~ = ve~ vm~ \n")
if mode == "W":
  f.write(f"generate p p > ax W+ , (W+ > l+ vl) , (ax > g g) \n")
  f.write(f"add process p p > ax W-, (W- > l- vl~), (ax > g g) \n")
if mode == "Z":
  f.write(f"generate p p > ax Z , Z > l+ l- , ax > g g \n")
f.write(f"output Script_{slug} \n")
f.write(f"launch Script_{slug} \n")
f.write(f"shower=Pythia8 \n") #Add Pythia
f.write(f"0 \n") #Launch the computation
f.write(f"set nevents = {nevent} \n" ) # change the number of event
f.write(f"set Ma {mass_Alp} \n" ) 
f.write(f"set CGtil 1 \n" )
f.write(f"set CWtil 0.01 \n" )
f.write(f"set CaPhi 0 \n" )
f.write(f"set CBtil -0.3049704673640523 \n" )
f.write(f"set decay 9000005 auto \n" )
f.write(f"set alpsfact 1 \n" )
f.write(f"set time_of_flight 0 \n" ) # Set the time of flight of the particle.
f.write(f"set event_norm = sum \n" )
f.write(f"set ptj = 0 \n" ) # Disable the cuts.
f.write(f"set ptb = 0 \n" )
f.write(f"set pta = 0 \n" )
f.write(f"set ptl = 0 \n" )
f.write(f"set etaj = -1 \n" )
f.write(f"set etab = -1 \n" )
f.write(f"set etaa = -1 \n" )
f.write(f"set etal = -1 \n" )
f.write(f"set drjj = 0 \n" )
f.write(f"set drbb = 0 \n" )
f.write(f"set drll = 0 \n" )
f.write(f"set draa = 0 \n" )
f.write(f"set drbj = 0 \n" )
f.write(f"set draj = 0 \n" )
f.write(f"set drjl = 0 \n" )
f.write(f"set drab = 0 \n" )
f.write(f"set drbl = 0 \n" )
f.write(f"set dral = 0 \n" )
f.write(f"0 \n") #launch the generation
f.write(f"exit")
