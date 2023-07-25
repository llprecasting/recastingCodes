import os

masses_Zp = [300,150,110,80,40,30]
masses_S = [200,100,80,50,5,1]

for i in range(len(masses_Zp)):
        f = open(f"Script_Asymetric_mZp{masses_Zp[i]}_mS{masses_S[i]}.txt", 'w')
        f.write(f"import model ./HAHM_variableMW_v3_UFO \n")
        f.write(f"define q = u c d s u~ c~ d~ s~ b b~ t t~ \n")
        f.write(f"define f = u c d s u~ c~ d~ s~ b b~ e+ e- mu+ mu- ta+ ta- t t~ \n")
        f.write(f"generate q q > Zp h2 HIG=1 HIW=1 QED=2 QCD=0,  (Zp > f f), (h2 > f f) \n")
        f.write(f"output Script_Asymetric_mZp{masses_Zp[i]}_mS{masses_S[i]} \n")
        f.write(f"launch Script_Asymetric_mZp{masses_Zp[i]}_mS{masses_S[i]} \n")
        f.write(f"1 \n")
        f.write(f"set pdlabel = lhapdf \n" )
        f.write(f"set lhaid = 315000 \n" )
        f.write(f"set mhsinput {masses_S[i]} \n") # dark Higgs
        f.write(f"set mZDinput {masses_Zp[i]} \n") # Dark Z
        f.write(f"set epsilon 1e-10 \n") # ATLAS : 1e-10
        f.write(f"set kap 1e-4 \n") # ATLAS : 1e-4
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
        f.write(f"set lhe_version = 3.0\n" )
        f.write(f"set cut_decays = F\n" )
        f.write(f"exit")


for mass,mass_s in zip(masses_Zp,masses_S):
    os.system(f"./bin/mg5_aMC -f Script_Asymetric_mZp{mass}_mS{mass_s}.txt")

