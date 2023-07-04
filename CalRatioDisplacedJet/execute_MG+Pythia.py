# This code write script that MG and Pythia will read to generate events
import os

masses_phi = [1000,600,400,200,125,60]
masses_S = [275,150,100,50,55,5]


for i in range(len(masses_phi)):
    if masses_phi[i] <= 125:
        f = open(f"script_mH{masses_phi[i]}_mS{masses_S[i]}.txt", 'w') # Creation of the Script to follow for MG.
        f.write(f"import model ./HAHM_variableMW_v3_UFO \n") # Import the model.
        f.write(f"define f = u c d s u~ c~ d~ s~ b b~ e+ e- mu+ mu- ta+ ta- t t~ \n") # Define a fermion.
        f.write(f"generate g g > h HIG=1 HIW=0 QED=0 QCD=0, (h > h2 h2, h2 > f f) \n") # Generate the process studied.
        f.write(f"output Script_mH{masses_phi[i]}_mS{masses_S[i]} \n")
        f.write(f"launch Script_mH{masses_phi[i]}_mS{masses_S[i]} \n")
        f.write(f"1 \n") #Add Pythia
        f.write(f"set pdlabel = lhapdf \n" ) # Set the Parton Distribution Function.
        f.write(f"set lhaid = 315000 \n" )
        f.write(f"set mhsinput {masses_S[i]} \n") # Set a mass for the LLP.
        f.write(f"set mhinput {masses_phi[i]} \n") # Set a mass for the Heavy Neutral Boson.
        f.write(f"set epsilon 1e-10 \n") # Set the couplings for epsilon.
        f.write(f"set kap 1e-4 \n") # Set the couplings for kappa.
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
        f.write(f"set cut_decays = F\n" )
        f.write(f"set use_syst = T \n" ) # Allows the code to keep track of the various parameter needed for the computation of the systematics
        f.write(f"set sys_scalefact = 1 0.5 2 \n" ) # Factorization/renormalization scale factor
        f.write(f"set sys_pdf = NNPDF31_lo_as_0118 \n" )
        f.write(f"set wzp = Auto\n"  ) # Set the calculation of the Width.
        f.write(f"set wh = Auto\n" )
        f.write(f"set wt = Auto\n" )
        f.write(f"set whs = Auto\n" )
        f.write(f"set lhe_version = 3.0\n" )
        f.write(f"exit")

    else:
        f = open(f"script_mH{masses_phi[i]}_mS{masses_S[i]}.txt", 'w')
        f.write(f"import model ./HAHM_variableMW_v3_UFO \n")
        f.write(f"define f = u c d s u~ c~ d~ s~ b b~ e+ e- mu+ mu- ta+ ta- t t~ \n")
        f.write(f"generate g g > h HIG=1 HIW=0 QED=0 QCD=0, (h > h2 h2, h2 > f f) \n")
        f.write(f"output Script_mH{masses_phi[i]}_mS{masses_S[i]} \n")
        f.write(f"launch Script_mH{masses_phi[i]}_mS{masses_S[i]} \n")
        f.write(f"1 \n")
        f.write(f"set pdlabel = lhapdf \n" )
        f.write(f"set lhaid = 315000 \n" )
        f.write(f"set mhsinput {masses_S[i]} \n")
        f.write(f"set mhinput {masses_phi[i]} \n")
        f.write(f"set epsilon 1e-10 \n")
        f.write(f"set kap 1e-4 \n")
        f.write(f"set time_of_flight 0 \n" )
        f.write(f"set event_norm = sum \n" )
        f.write(f"set ptj = 0 \n" )
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
        f.write(f"set use_syst = T \n" )
        f.write(f"set sys_scalefact = 1 0.5 2 \n" )
        f.write(f"set sys_pdf = NNPDF31_lo_as_0118 \n" )
        f.write(f"set wzp = 5\n" )
        f.write(f"set wh = 5\n" )
        f.write(f"set wt = Auto\n" )
        f.write(f"set whs = 5\n" )
        f.write(f"exit")

for mass_phi,mass_s in zip(masses_phi,masses_S):
    os.system(f"./bin/mg5_aMC -f script_mH{mass_phi}_mS{mass_s}.txt") # Launch the Script in MG.

