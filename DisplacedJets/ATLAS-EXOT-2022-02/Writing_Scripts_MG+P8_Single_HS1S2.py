# This code writes script that MG and Pythia will read to generate events
import os
import sys


mode = sys.argv[1]
mass_pair = sys.argv[2]
mass_Phi = int(mass_pair.split("_")[0])
mass_S1 = int(mass_pair.split("_")[1])
mass_S2 = int(mass_pair.split("_")[2])



nevent = int(sys.argv[3])
slug = f'{mode}H{mass_Phi}_S1{mass_S1}_S2{mass_S2}_{mode}' 
InDir = sys.argv[4]
OutDir = sys.argv[5]
Lambda ="100"
k="1"
#ct1 = sys.argv[9]
#ct2= sys.argv[10]

f = open(f"{OutDir}/script_{slug}.txt", 'w') # Creation of the Script to follow for MG.
f.write(f"import model sm \n")
f.write(f"import model /users/divers/atlas/millot/home2/MG5_aMC_v3_4_2/HAHM_MG5model_v3/HAHM_asymmetric_UFO/ \n") # Import the model.
f.write(f"define f = u c d s u~ c~ d~ s~ b b~ e+ e- mu+ mu- ta+ ta- t t~ \n") # Define a fermion.
f.write(f"generate g g > h NP=99 HIG=99 QED=99 QCD=99, (h > 35 36, 35 > f f, 36 > f f) \n") # Generate the process.
f.write(f"output Script_{slug} \n")
f.write(f"launch Script_{slug} \n")
f.write(f"shower=Pythia8 \n") #Add Pythia
f.write(f"0 \n") #Launch the computation
f.write(f"set nevents = {nevent} \n" ) # change the number of event
f.write(f"set ms1 {mass_S1} \n") # Set a mass for the LLP.
f.write(f"set ms2 {mass_S2} \n") # Set a mass for the LLP.
f.write(f"set mh {mass_Phi} \n") # Set a mass for the Heavy Neutral Boson.
#f.write(f"set ct1 1.0 \n") # Set a mass for the LLP.
#f.write(f"set ct2 1.0 \n") # Set a mass for the LLP.
f.write(f"set ct1 1000 \n") # Set a mass for the LLP.
f.write(f"set ct2 0.001 \n") # Set a mass for the LLP.
f.write(f"set Lambda {Lambda} \n")
f.write(f"set k {k}\n") # Set the couplings for kappa.
f.write(f"set cl 1 \n")
f.write(f"set time_of_flight 0 \n" ) # Set the time of flight of the particle.
f.write(f"set event_norm = sum \n" )
f.write(f"set ptj = 0 \n" ) # Disable the cuts.
f.write(f"set ptb = 0 \n" )
f.write(f"set pta = 0 \n" )
f.write(f"set ptl = 0 \n" )
f.write(f"set iseed = 12346 \n" )
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
f.write(f"set use_syst = T\n") 
f.write(f"set wh = Auto\n")
f.write(f"set Ws1 = Auto\n")
f.write(f"set Ws2 = Auto\n")

if mass_Phi <= 125: 
#if False: 
    f.write(f"set cut_decays = F\n" )
    f.write(f"set wh = 7.290534e-03\n" )
    f.write(f"set wt = Auto\n" )
#    f.write(f"set whs = Auto\n" )
else:
 #   f.write(f"set use_syst = T \n" )
    f.write(f"set wzp = 5\n" )
    f.write(f"set wh = 5\n" )
    f.write(f"set wt = Auto\n" )
 #   f.write(f"set whs = 5\n" )

f.write(f"0 \n") #launch the generation
f.write(f"exit")
