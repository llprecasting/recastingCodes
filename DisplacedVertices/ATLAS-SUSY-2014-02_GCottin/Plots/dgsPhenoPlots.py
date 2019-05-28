#Plotting script for displaced dgsPheno project
#Written by Giovanna Cottin
import os,sys
import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.random import normal
import matplotlib.font_manager as font_manager
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter 
from pylab import *
import matplotlib.gridspec as gridspec
###########
#Data files
###########
effVsLifetime_GGM          = "DataFiles/effVsCtau_GGM_8TeV_100k.dat"
effVsLifetime_GGM2         = "DataFiles/effVsCtau_GGM2_8TeV_100k.dat"
effVsLifetime_RPV          = "DataFiles/effVsCtau_RPV_8TeV_100k.dat"
mDVnTrk                    = "DataFiles/mDVnTrk_DGS_8TeV_100k.dat"
mDVnTrk2                   = "DataFiles/mDVnTrk_DGSA170_8TeV_100k.dat"
effVsA1Mass                = "DataFiles/effVsA1Mass_DGS_8TeV_100k.dat"
effVsLifetime_DGS   	     = "DataFiles/effVsCtau_DGS_8TeV_500k.dat"
effVsLifetime_DGS_13TeV    = "DataFiles/effVsCtau_DGS_13TeV_500k.dat"

effVsLifetime_DGS_5510     = "DataFiles/effVsCtau_DGS_5510_8TeV_500k.dat"
effVsLifetime_DGS_555      = "DataFiles/effVsCtau_DGS_555_8TeV_500k.dat"
effVsLifetime_DGS_525      = "DataFiles/effVsCtau_DGS_525_8TeV_500k.dat"
effVsLifetime_DGS_5256jtp  = "DataFiles/effVsCtau_DGS_525prompt_8TeV_500k.dat"

effVsLifetime_DGS_5256jtp_13TeV  = "DataFiles/effVsCtau_DGS_525prompt_DVMass_13TeV_100k_NEW.dat"

effVsLifetime_DGS_5256jtp_13TeV_Jets    = "DataFiles/effVsCtau_DGS_525prompt_Jets_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_Meff    = "DataFiles/effVsCtau_DGS_525prompt_Meff_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_MET     = "DataFiles/effVsCtau_DGS_525prompt_MET_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_Dphi    = "DataFiles/effVsCtau_DGS_525prompt_dphi_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_Ratio   = "DataFiles/effVsCtau_DGS_525prompt_ratio_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_DV      = "DataFiles/effVsCtau_DGS_525prompt_DV_13TeV_100k_NEW.dat"
effVsLifetime_DGS_5256jtp_13TeV_DVNtrk  = "DataFiles/effVsCtau_DGS_525prompt_DVNtrk_13TeV_100k_NEW.dat"


##########################
#Data for ATLAS Benchmarks
###############################################
# ATLAS GGM mchi = 400 GeV ; m_gluino = 1.1 TeV
###############################################
nominalX = np.array([1,1.43845,2.06914,2.97635,4.28133,6.15848,8.85867,12.7427,18.3298,26.3665,37.9269,54.5559,78.476,112.884,162.378,233.572,335.982,483.293,695.193,1000]) 
nominalY=np.array([0.000384428,0.00337376,0.0147529,0.0386017,0.071398,0.106704,0.140067,0.167443,0.185385,0.192887,0.190499,0.178937,0.15931,0.133699,0.105173,0.0772114,0.0528409,0.0338494,0.0204674,0.0117909])
Err=np.array([[1,0.000624983],[1.43845,0.00500215],[2.06914,0.0206813],[2.97635,0.0521283],[4.28133,0.0939334],[6.15848,0.138166],[8.85867,0.180561],[12.7427,0.216596],[18.3298,0.24089],[26.3665,0.251396],[37.9269,0.248825],[54.5559,0.234204],[78.476,0.208935],[112.884,0.175658],[162.378,0.138367],[233.572,0.102017],[335.982,0.0706353],[483.293,0.0469472],[695.193,0.0308608],[1000,0.0204484],[1000,0.00768328],[695.193,0.0136653],[483.293,0.0229349],[335.982, 0.0360273],[233.572,0.0527525],[162.378,0.0719785],[112.884,0.0917403],[78.476,0.109685],[54.5559,0.12367],[37.9269,0.132173],[26.3665,0.134379],[18.3298,0.12988],[12.7427,0.118291],[8.85867,0.0995729],[6.15848,0.075243],[4.28133, 0.0488626],[2.97635,0.0250751],[2.06914,0.00882454],[1.43845,0.00174537],[1,0.000143872]])
##############################################
# ATLAS GGM2 mchi = 1 TeV ; m_gluino = 1.1 TeV
##############################################
nominalGGM2X=np.array([1,1.43845,2.06914,2.97635,4.28133,6.15848,8.85867,12.7427,18.3298,26.3665,37.9269,54.5559,78.476,112.884,162.378,233.572,335.982,483.293,695.193,1000])                            
nominalGGM2Y=np.array([0.00027905,0.00158573,0.00444047,0.00837391,0.0138648,0.0223146,0.0337031,0.0460647,0.0568487,0.0639989,0.0662822,0.0638444,0.0580141,0.0503075,0.0418354,0.0333067,0.0253138,0.018363,0.0126699,0.00825567])
errGGM2=np.array([[1,0.00057638],[1.43845,0.00325747],[2.06914,0.00890491],[2.97635,0.0155874],[4.28133,0.0226468],[6.15848,0.0327787],[8.85867,0.0477734],[12.7427,0.0652375],[18.3298,0.0809187],[26.3665,0.091208],[37.9269,0.0942126],[54.5559,0.0903937],[78.476,0.0818278],[112.884,0.0707009],[162.378,0.0585797],[233.572,0.0465041],[335.982,0.0353082],[483.293,0.0256887],[695.193,0.0180194],[1000,0.0123358],[1000,0.00476993],[695.193,0.0075378],[483.293,0.0110833],[335.982,0.0153243],[233.572,0.02011],[162.378,0.0250912],[112.884,0.0299141],[78.476,0.0342003],[54.5559,0.0372952],[37.9269,0.0383517],[26.3665,0.0367897],[18.3298,0.0327787],[12.7427,0.026892],[8.85867,0.0196329],[6.15848,0.0118505],[4.28133,0.00508285],[2.97635,0.00116043],[2.06914,0],[1.43845,0],[1, 0]])
##############################################
#ATLAS RPV mchi = 500 GeV ; m_squark = 700 GeV
##############################################
rpvNominalX=np.array([1,1.43845,2.06914,2.97635,4.28133,6.15848,8.85867,12.7427,18.3298,26.3665,37.9269,54.5559,78.476,112.884,162.378,233.572,335.982,483.293,695.193,1000])           
rpvNominalY=np.array([0.00408635,0.020695,0.0616625,0.129378,0.214515,0.302024,0.378343,0.435519,0.471063,0.487041,0.487238,0.473776,0.447386,0.409167,0.361299,0.30715,0.250353,0.193498,0.139644,0.0932954])
rpvErr=np.array([[1,0.00549255],[1.43845,0.0262213],[2.06914,0.0738698],[2.97635,0.147597],[4.28133,0.238275],[6.15848,0.331988],[8.85867,0.414163],[12.7427,0.475199],[18.3298,0.512501],[26.3665,0.528913],[37.9269,0.528625],[54.5559,0.51361],[78.476,0.484453],[112.884,0.442379],[162.378,0.389873],[233.572,0.330769],[335.982,0.26955],[483.293,0.210944],[695.193,0.161097],[1000,0.12243],[1000,0.0852041],[695.193,0.128627],[483.293,0.178976],[335.982,0.231584],[233.572,0.283561],[162.378,0.332725],[112.884,0.375956],[78.476,0.41032],[54.5559,0.433941],[37.9269,0.44585],[26.3665,0.445168],[18.3298,0.429626],[12.7427,0.39584],[8.85867,0.342523],[6.15848,0.272059],[4.28133,0.190755],[2.97635,0.111159],[2.06914,0.0494551],[1.43845,0.0151687],[1,0.00268015]])
###########################
#Extracting info from files
###########################
effVsLifetime_GGM_data                = {0:[],1:[],2:[]}
effVsLifetime_GGM2_data               = {0:[],1:[],2:[]}
effVsLifetime_RPV_data                = {0:[],1:[],2:[]}
effVsLifetime_DGS_data                = {0:[],1:[],2:[]}
effVsLifetime_DGS_data_13TeV          = {0:[],1:[],2:[]}
effVsLifetime_DGS_5510_data           = {0:[],1:[],2:[]}
effVsLifetime_DGS_555_data            = {0:[],1:[],2:[]}
effVsLifetime_DGS_525_data            = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data        = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV  = {0:[],1:[],2:[]}


effVsLifetime_DGS_5256jtp_data_13TeV_Jets   = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_Meff   = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_MET    = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_Dphi   = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_Ratio  = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_DV     = {0:[],1:[],2:[]}
effVsLifetime_DGS_5256jtp_data_13TeV_DVNtrk = {0:[],1:[],2:[]}


mDVnTrk_data                          = {0:[],1:[]}
mDVnTrk2_data                         = {0:[],1:[]}
effVsA1Mass_data                      = {0:[],1:[],2:[]}

for line in open(effVsLifetime_GGM,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_GGM_data[i].append(float(C))
for line in open(effVsLifetime_GGM2,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_GGM2_data[i].append(float(C))
for line in open(effVsLifetime_RPV,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_RPV_data[i].append(float(C))
for line in open(mDVnTrk,"r"):
  for i,C in enumerate(line.split()): mDVnTrk_data[i].append(float(C))
for line in open(mDVnTrk2,"r"):
  for i,C in enumerate(line.split()): mDVnTrk2_data[i].append(float(C))
for line in open(effVsA1Mass,"r"):
  for i,C in enumerate(line.split()): effVsA1Mass_data[i].append(float(C))
for line in open(effVsLifetime_DGS,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_data[i].append(float(C))
for line in open(effVsLifetime_DGS_5510,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5510_data[i].append(float(C))
for line in open(effVsLifetime_DGS_555,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_555_data[i].append(float(C))
for line in open(effVsLifetime_DGS_525,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_525_data[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data[i].append(float(C))
for line in open(effVsLifetime_DGS_13TeV,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_data_13TeV[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV[i].append(float(C))


for line in open(effVsLifetime_DGS_5256jtp_13TeV_Jets,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_Jets[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_Meff,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_Meff[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_MET,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_MET[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_Dphi,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_Dphi[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_Ratio,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_Ratio[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_DV,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_DV[i].append(float(C))
for line in open(effVsLifetime_DGS_5256jtp_13TeV_DVNtrk,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_DGS_5256jtp_data_13TeV_DVNtrk[i].append(float(C))


ctau_ATLASBM         = np.array(effVsLifetime_GGM_data[0])
effGGM               = np.array(effVsLifetime_GGM_data[2])
effGGM2              = np.array(effVsLifetime_GGM2_data[2])
effRPV               = np.array(effVsLifetime_RPV_data[2])
effDGS               = np.array(effVsLifetime_DGS_data[2])
effDGS_5510          = np.array(effVsLifetime_DGS_5510_data[2])
effDGS_555           = np.array(effVsLifetime_DGS_555_data[2])
effDGS_525           = np.array(effVsLifetime_DGS_525_data[2])
effDGS_5256jtp       = np.array(effVsLifetime_DGS_5256jtp_data[2])
effDGS_13TeV         = np.array(effVsLifetime_DGS_data_13TeV[2])
effDGS_5256jtp_13TeV = np.array(effVsLifetime_DGS_5256jtp_data_13TeV[2])

effDGS_5256jtp_13TeV_Jets = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_Jets[2])
effDGS_5256jtp_13TeV_Meff = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_Meff[2])
effDGS_5256jtp_13TeV_MET = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_MET[2])
effDGS_5256jtp_13TeV_Dphi = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_Dphi[2])
effDGS_5256jtp_13TeV_Ratio = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_Ratio[2])
effDGS_5256jtp_13TeV_DV = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_DV[2])
effDGS_5256jtp_13TeV_DVNtrk = np.array(effVsLifetime_DGS_5256jtp_data_13TeV_DVNtrk[2])


mass                 = np.array(mDVnTrk_data[0])
mass2                = np.array(mDVnTrk2_data[0])
tracks               = np.array(mDVnTrk_data[1])
tracks2              = np.array(mDVnTrk2_data[1])
A1mass               = np.array(effVsA1Mass_data[0])
effA1                = np.array(effVsA1Mass_data[2])

###########################
#Computing stat error bands
###########################
nPassGGM     = np.array(effVsLifetime_GGM_data[1])
nPassGGM2    = np.array(effVsLifetime_GGM2_data[1])
nPassRPV     = np.array(effVsLifetime_RPV_data[1])
nPassDGS            = np.array(effVsLifetime_DGS_data[1])
nPassDGS_5510       = np.array(effVsLifetime_DGS_5510_data[1])
nPassDGS_555        = np.array(effVsLifetime_DGS_555_data[1])
nPassDGS_525        = np.array(effVsLifetime_DGS_525_data[1])
nPassDGS_5256jtp    = np.array(effVsLifetime_DGS_5256jtp_data[1])

nPassDGS_13TeV         = np.array(effVsLifetime_DGS_data_13TeV[1])
nPassDGS_5256jtp_13TeV = np.array(effVsLifetime_DGS_5256jtp_data_13TeV[1])
nPasseffA1   = np.array(effVsA1Mass_data[1])

deltaEffGGM  = 1.0/np.sqrt(nPassGGM)
deltaEffGGM2 = 1.0/np.sqrt(nPassGGM2)
deltaEffRPV  = 1.0/np.sqrt(nPassRPV)
deltaEffA1   = 1.0/np.sqrt(nPasseffA1) 
deltaEffDGS = 1.0/np.sqrt(nPassDGS) 
deltaEffDGS_13TeV = 1.0/np.sqrt(nPassDGS_13TeV) 

deltaEffDGS_5510 = 1.0/np.sqrt(nPassDGS_5510) 
deltaEffDGS_555 = 1.0/np.sqrt(nPassDGS_555) 
deltaEffDGS_525 = 1.0/np.sqrt(nPassDGS_525) 
deltaEffDGS_5256jtp = 1.0/np.sqrt(nPassDGS_5256jtp) 
deltaEffDGS_5256jtp_13TeV = 1.0/np.sqrt(nPassDGS_5256jtp_13TeV) 

# Stat error within 10%, nPass at least 100
# deltaEff=0.1

#Errors, how away we are from ATLAS
errGGMUp= Err[...,1][0:20][::-1]#from big to low ctau we saved our eff tables
errGGMDown= Err[...,1][20::]
ATLASGGM_error1=errGGMUp-nominalY[::-1]
ATLASGGM_error2=nominalY[::-1]-errGGMDown
errGGM2Up=errGGM2[...,1][0:20][::-1]
errGGM2Down= errGGM2[...,1][20::]
ATLASGGM2_error1=errGGM2Up-nominalGGM2Y[::-1]
ATLASGGM2_error2=nominalGGM2Y[::-1]-errGGM2Down
errRPVUp=rpvErr[...,1][0:20][::-1]
errRPVDown= rpvErr[...,1][20::]
ATLASRPV_error1=errRPVUp-rpvNominalY[::-1]
ATLASRPV_error2=rpvNominalY[::-1]-errRPVDown
#errors are almost symmetric, except slightly not on edges, taking simple average
ATLASGGM_error  = (np.abs(ATLASGGM_error1+ATLASGGM_error2))/2
ATLASGGM2_error = (np.abs(ATLASGGM2_error1+ATLASGGM2_error2))/2
ATLASRPV_error  = (np.abs(ATLASRPV_error1+ATLASRPV_error2))/2
#how far away are we from ATLAS results?
diffGGM  =(effGGM-nominalY[::-1])/ATLASGGM_error
diffGGM2 =(effGGM2-nominalGGM2Y[::-1])/ATLASGGM2_error
diffRPV  =(effRPV-rpvNominalY[::-1])/ATLASRPV_error

print "max diffGGM : ", np.max(diffGGM)
print "max diffGGM2 : ", np.max(diffGGM2)
print "max diffRPV : ", np.max(diffRPV)

###############
#Plotting Style
###############
rcParams['legend.loc'] = 'best'
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 20,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'legend.fontsize': 12
          # 'legend.fontsize': 20
          }
plt.rcParams.update(params) 

lp=0.27
rp=0.75
tp=0.94
bp=0.14
#modify bottom y left
#for n subplots is /(n+(n-1)wp  
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
# fig1=plt.figure(1,figsize=(2*(3.37*2), 2*(3.37*2*plotaspect/1.62)))
fig1=plt.figure(1,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
# I create groups that share attributes, in this case, I will have two groups with only 1 subplot each
gs1 = gridspec.GridSpec(1, 1) # here I define the amount of subplots I want in this subgroup
gs2 = gridspec.GridSpec(1, 1)
# The first plot, gs1, is gonna be from 0.95 and 0.3 in y (that full range is 0 to 1)
gs1.update(top=tp, bottom=0.3) # subplot 1 attributes
gs2.update(top=0.3,bottom=bp)  # subplot 2 attributes
#plt.subplot(gs1[0, 0]) # now we plot
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=25)
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.xlim([1., 1000.])
plt.ylim([0.0, 0.45])
plt.xscale('log')
plt.plot(nominalX, nominalY, '-b', alpha=0.5, color='b',linestyle='-', label="ATLAS DV+jets GGM")
plt.plot(Err[...,0][0:20], Err[...,1][0:20], '-b', alpha=0.5)
plt.plot(Err[...,0][20::], Err[...,1][20::], '-b', alpha=0.5)
plt.fill_between(nominalX, Err[...,1][0:20],Err[...,1][20::][::-1], color='b', edgecolor='b', alpha=0.05)
plt.plot(ctau_ATLASBM, effGGM, color='b',linestyle='--',lw=2,  label=" Our Simulation 8 TeV")
plt.fill_between(ctau_ATLASBM,effGGM-deltaEffGGM*effGGM,effGGM+deltaEffGGM*effGGM,color='b',alpha=0.3)
lg = plt.legend()
lg.draw_frame(False)
plt.ylabel('event efficiency',fontsize=30)
plt.gca().axes.xaxis.set_ticklabels([])
#hide the 0
plt.gcf().canvas.draw()
A=plt.gca().axes.yaxis.get_ticklabels()
B=[y.get_text() for y in A]
B[0]=u""
plt.gca().axes.yaxis.set_ticklabels(B)
ax2=plt.subplot(gs2[0, 0]) # now we plot
ax2.tick_params(axis='x', labelsize=25)
ax2.tick_params(axis='y', labelsize=10)
plt.plot(ctau_ATLASBM,diffGGM,color='b',lw=2,alpha=0.3)
plt.axhline(y=0,color='k',linestyle='--',lw=2,alpha=0.3)
plt.xlim([1., 1000.])
plt.ylim(-5.0,5.0)
plt.xscale('log')
plt.ylabel('$\\frac{\\textrm{eff} - \\textrm{ATLAS eff}}{\\textrm{ATLAS error}}$',fontsize=15)
plt.xlabel('$c\\tau$ [mm]',fontsize=30)

plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
# fig2=plt.figure(2,figsize=(2*(3.37*2), 2*(3.37*2*plotaspect/1.62)))
fig2=plt.figure(2,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3)
gs2.update(top=0.3,bottom=bp)
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=25)
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.xlim([1., 1000.])
plt.ylim([0.0, 0.45])
plt.xscale('log')
plt.plot(nominalGGM2X, nominalGGM2Y, '-m', alpha=0.5, color='m',linestyle='-', label="ATLAS DV+jets GGM2")
plt.plot(errGGM2[...,0][0:20], errGGM2[...,1][0:20], '-m', alpha=0.5)
plt.plot(errGGM2[...,0][20::], errGGM2[...,1][20::], '-m', alpha=0.5)
plt.fill_between(nominalGGM2X, errGGM2[...,1][0:20],errGGM2[...,1][20::][::-1], color='m', edgecolor='m', alpha=0.05)
plt.plot(ctau_ATLASBM, effGGM2, color='m',linestyle='--',lw=2,  label=" Our Simulation 8 TeV")
plt.fill_between(ctau_ATLASBM,effGGM2-deltaEffGGM2*effGGM2,effGGM2+deltaEffGGM2*effGGM2,color='m',alpha=0.3)
lg = plt.legend()
lg.draw_frame(False)
plt.ylabel('event efficiency',fontsize=30)
plt.gca().axes.xaxis.set_ticklabels([])
plt.gcf().canvas.draw()
A=plt.gca().axes.yaxis.get_ticklabels()
B=[y.get_text() for y in A]
B[0]=u""
plt.gca().axes.yaxis.set_ticklabels(B)
ax2=plt.subplot(gs2[0, 0]) # now we plot
ax2.tick_params(axis='x', labelsize=25)
ax2.tick_params(axis='y', labelsize=10)
plt.plot(ctau_ATLASBM,diffGGM2,color='m',lw=2,alpha=0.3)
plt.axhline(y=0,color='k',linestyle='--',lw=2,alpha=0.3)
plt.xlim([1., 1000.])
plt.ylim(-5.0,5.0)
plt.xscale('log')
plt.ylabel('$\\frac{\\textrm{eff} - \\textrm{ATLAS eff}}{\\textrm{ATLAS error}}$',fontsize=15)
plt.xlabel('$c\\tau$ [mm]',fontsize=30)

plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig3=plt.figure(3,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3) 
gs2.update(top=0.3,bottom=bp)  
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=25)
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.xlim([1., 1000.])
plt.ylim([0.0, 1.0])
plt.xscale('log')
plt.plot(rpvNominalX, rpvNominalY, '-g', alpha=0.5, color='g',linestyle='-', label="ATLAS DV+jets RPV")
plt.plot(rpvErr[...,0][0:20], rpvErr[...,1][0:20], '-g', alpha=0.5)
plt.plot(rpvErr[...,0][20::], rpvErr[...,1][20::], '-g', alpha=0.5)
plt.fill_between(rpvNominalX, rpvErr[...,1][0:20],rpvErr[...,1][20::][::-1], color='g', edgecolor='g', alpha=0.05)
plt.plot(ctau_ATLASBM, effRPV, color='g',linestyle='--',lw=2,  label=" Our Simulation 8 TeV")
plt.fill_between(ctau_ATLASBM,effRPV-deltaEffRPV*effRPV,effRPV+deltaEffRPV*effRPV,color='g',alpha=0.3)
lg = plt.legend()
lg.draw_frame(False)
plt.ylabel('event efficiency',fontsize=30)
plt.gca().axes.xaxis.set_ticklabels([])
plt.gcf().canvas.draw()
A=plt.gca().axes.yaxis.get_ticklabels()
B=[y.get_text() for y in A]
B[0]=u""
plt.gca().axes.yaxis.set_ticklabels(B)
ax2=plt.subplot(gs2[0, 0]) # now we plot
ax2.tick_params(axis='x', labelsize=25)
ax2.tick_params(axis='y', labelsize=10)
plt.plot(ctau_ATLASBM,diffRPV,color='g',lw=2,alpha=0.3)
plt.axhline(y=0,color='k',linestyle='--',lw=2,alpha=0.3)
plt.xlim([1., 1000.])
plt.ylim(-5.0,5.0)
plt.xscale('log')
plt.ylabel('$\\frac{\\textrm{eff} - \\textrm{ATLAS eff}}{\\textrm{ATLAS error}}$',fontsize=15)
plt.xlabel('$c\\tau$ [mm]',fontsize=30)

####################################
#Fraction of events in signal region
####################################
x0=5
x1=10
y0=10.0
y1=20.0

#[0] since np.where returns more than one thing, but [0] is the element with the indices with the conditions
signalEvents1=len(np.where((tracks>x0) & (tracks<x1) & (mass>y0) & (mass<y1))[0])
signalEvents2=len(np.where((tracks2>x0) & (tracks2<x1) & (mass2>y0) & (mass2<y1))[0])
# print"signalEvents1 = ",signalEvents1
# print"signalEvents2 = ",signalEvents2
fraction1=signalEvents1*100.0/len(mass)
fraction2=signalEvents2*100.0/len(mass2)
print "fraction1 = ",fraction1
print "fraction2 = ",fraction2

lp=0.1
rp=0.97
tp=0.90
bp=0.18
wp=0.0

#for n subplots is /(n+(n-1)wp  
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp))/(2+wp) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
#fig4=plt.figure(4,figsize=(2*(3.37*2), 2*(3.37*2*plotaspect/1.62)))
fig4=plt.figure(4,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp,wspace=wp)
ax4 = fig4.add_subplot(121)
# histogram2d binning gives: locx=np.linspace(np..min(x),np.max(x),endpoint=True)
h1 , locx, locy = np.histogram2d(tracks,mass,bins=[6.,20.])#,range=np.array([(2.,10.), (0.1, 20.)]))
h1=h1.T# histogram2d shows this inverted
plt.imshow(h1,interpolation='none',origin='lower',extent=[locx[0],locx[-1],locy[0],locy[-1]],aspect='auto',cmap='binary',norm=LogNorm())
plt.ylim(0)
plt.ylim([0,20])
plt.xlim([2,9.8])
plt.plot([5,10],[10,10],color='m',lw=2)
plt.plot([5,5],[10,20],color='m',lw=2)
ax4.text(7.3,13,"Signal Region",fontsize=15,color='m')
ax4.text(6.8,11,"efficiency of $0.01\\%$",fontsize=15,color='k')
ax4.text(6.8,3,"$c\\tau_{\\tilde{\chi}^{0}_{1}}=64$ mm ",fontsize=20)
ax4.text(6.8,1,"$m_{a_{1}}=26$ GeV ",fontsize=20)
plt.ylabel('DV invariant mass [GeV]',fontsize=30)
plt.xlabel('DV track multiplicity',fontsize=30)
ax4.tick_params(axis='x', labelsize=25)
ax4.tick_params(axis='y', labelsize=25)
ax5 = fig4.add_subplot(122)
#In order to have same binning as fig4. locx and locy are the positions of the boarder of the binns 
h2 , locx2, locy2 = np.histogram2d(tracks2,mass2,bins=[locx,locy])
h2=h2.T# histogram2d shows this inverted
plt.imshow(h2,interpolation='none',origin='lower',extent=[locx2[0],locx2[-1],locy2[0],locy2[-1]],aspect='auto',cmap='binary',norm=LogNorm())
plt.ylim(0)
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.ylim([0,20])
plt.xlim([2,10])
plt.plot([5,10],[10,10],color='m',lw=2)
plt.plot([5,5],[10,20],color='m',lw=2)
ax5.text(7.5,13,"Signal Region",fontsize=15,color='m')
ax5.text(7.0,11,"efficiency of $0.15\\%$",fontsize=15,color='k')
ax5.text(7.0,3,"$c\\tau_{\\tilde{\chi}^{0}_{1}}=64$ mm ",fontsize=20)
ax5.text(7.0,1,"$m_{a_{1}}=70$ GeV ",fontsize=20)
plt.gca().axes.yaxis.set_ticklabels([])
plt.xlabel('DV track multiplicity',fontsize=30)
ax5.tick_params(axis='x', labelsize=25)
ax5.tick_params(axis='y', labelsize=25)

lp=0.27
rp=0.75
tp=0.94
bp=0.16
#modify bottom y left
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
#fig6=plt.figure(6,figsize=(2*(3.37*2), 2*(3.37*2*plotaspect/1.62)))
fig6=plt.figure(6,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
plt.xlim([20., 100.])
# plt.ylim([0.0, 0.02])
plt.yscale('log')
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.plot(A1mass, effA1, color='r',linestyle='-.',lw=2,alpha=0.5)
plt.fill_between(A1mass,effA1-deltaEffA1*effA1,effA1+deltaEffA1*effA1,color='r',alpha=0.3)
ax6 = fig6.add_subplot(111)
ax6.tick_params(axis='x', labelsize=25)
ax6.tick_params(axis='y', labelsize=25)
plt.ylabel('event efficiency',fontsize=30)
plt.xlabel('$m_{a_{1}}$ [GeV]',fontsize=30)

lp=0.27
rp=0.75
tp=0.94
bp=0.16
#modify bottom y left
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig7=plt.figure(7,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
plt.xlim([1., 1000.])
plt.yscale('log')
plt.xscale('log')
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.plot(ctau_ATLASBM, effDGS_525, color='b',linestyle='-',lw=2,alpha=0.5, label= "5 mm+ 2 trk  + 5 GeV inv. mass ")
plt.plot(ctau_ATLASBM, effDGS_555, color='c',linestyle='--',lw=2,alpha=0.5, label= "5 mm+ 5 trk  + 5 GeV inv. mass ")
plt.plot(ctau_ATLASBM, effDGS_5256jtp, color='m',linestyle=':',lw=2,alpha=0.5, label= "5 mm+ 2 trk  + 5 GeV inv. mass + prompt cuts ")
plt.plot(ctau_ATLASBM, effDGS_5510, color='k',linestyle='-.',lw=2,alpha=0.5, label= "5 mm+ 5 trk  + 10 GeV inv. mass ")
plt.plot(ctau_ATLASBM, effDGS, color='r',linestyle='--',lw=2,alpha=0.5, label= "1 mm+ 5 trk  + 10 GeV inv. mass (ATLAS) ")
plt.fill_between(ctau_ATLASBM,effDGS_525-deltaEffDGS_525*effDGS_525,effDGS_525+deltaEffDGS_525*effDGS_525,color='b',alpha=0.3)
plt.fill_between(ctau_ATLASBM,effDGS_555-deltaEffDGS_555*effDGS_555,effDGS_555+deltaEffDGS_555*effDGS_555,color='c',alpha=0.3)
plt.fill_between(ctau_ATLASBM,effDGS_5256jtp-deltaEffDGS_5256jtp*effDGS_5256jtp,effDGS_5256jtp+deltaEffDGS_5256jtp*effDGS_5256jtp,color='m',alpha=0.3)
plt.fill_between(ctau_ATLASBM,effDGS_5510-deltaEffDGS_5510*effDGS_5510,effDGS_5510+deltaEffDGS_5510*effDGS_5510,color='k',alpha=0.3)
plt.fill_between(ctau_ATLASBM,effDGS-deltaEffDGS*effDGS,effDGS+deltaEffDGS*effDGS,color='r',alpha=0.3)
ax6 = fig6.add_subplot(111)
ax6.tick_params(axis='x', labelsize=25)
ax6.tick_params(axis='y', labelsize=25)
plt.ylabel('event efficiency',fontsize=30)
plt.xlabel('$c\\tau$ [mm]',fontsize=30)
lg = plt.legend()
lg.draw_frame(False)

lp=0.27
rp=0.75
tp=0.94
bp=0.16
#modify bottom y left
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig8=plt.figure(8,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
plt.xlim([1., 1000.])
# plt.ylim([0.0, 0.01])
plt.yscale('log')
plt.xscale('log')
plt.title("Allanach, Badziak, Cottin, Desai, Hugonie, Ziegler (2016)", fontsize=15, loc='right')
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_Jets, color='b',linestyle='-',lw=2,alpha=0.5, label= "All jet cuts ")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_Meff, color='k',linestyle='-.',lw=2,alpha=0.5, label= " Meff")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_MET, color='c',linestyle='--',lw=2,alpha=0.5, label= " MET")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_Dphi, color='r',linestyle=':',lw=2,alpha=0.5, label= "Dphi ")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_Ratio, color='g',linestyle='-',lw=2,alpha=0.5, label= " MET/Meff")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_DV, color='y',linestyle='--',lw=2,alpha=0.5, label= "DVReco (in fidutial and material) ")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV_DVNtrk, color='b',linestyle='-.',lw=2, alpha=0.5,label= " DVNtrk ")
plt.plot(ctau_ATLASBM, effDGS_5256jtp_13TeV, color='m',linestyle=':',lw=2,alpha=0.5, label= "5 mm+ 2 trk  + 5 GeV inv. mass + prompt cuts  (or DVMass)")
# plt.fill_between(ctau_ATLASBM,effDGS_5256jtp-deltaEffDGS_5256jtp*effDGS_5256jtp,effDGS_5256jtp+deltaEffDGS_5256jtp*effDGS_5256jtp,color='m',alpha=0.3)

ax6 = fig6.add_subplot(111)
ax6.tick_params(axis='x', labelsize=25)
ax6.tick_params(axis='y', labelsize=25)
plt.ylabel('event efficiency',fontsize=30)
plt.xlabel('$c\\tau$ [mm]',fontsize=30)
lg = plt.legend()
lg.draw_frame(False)








# lp=0.27
# rp=0.75
# tp=0.94
# bp=0.14
# #modify bottom y left
# #for n subplots is /(n+(n-1)wp  
# plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
# # fig1=plt.figure(1,figsize=(2*(3.37*2), 2*(3.37*2*plotaspect/1.62)))
# fig8=plt.figure(8,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
# plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
# # I create groups that share attributes, in this case, I will have two groups with only 1 subplot each
# gs1 = gridspec.GridSpec(1, 1) # here I define the amount of subplots I want in this subgroup
# gs2 = gridspec.GridSpec(1, 1)
# # The first plot, gs1, is gonna be from 0.95 and 0.3 in y (that full range is 0 to 1)
# gs1.update(top=tp, bottom=0.3) # subplot 1 attributes
# gs2.update(top=0.3,bottom=bp)  # subplot 2 attributes
# #plt.subplot(gs1[0, 0]) # now we plot
# ax1=plt.subplot(gs1[0, 0])
# ax1.tick_params(axis='y', labelsize=25)
# plt.xlim([-8, 8.])
# plt.ylim([0.01, 2.2])
# plt.hist(d0Bs[abs(d0Bs)>0],color="black",range=[-8.0,8.0],histtype="step",lw=2,label=" $t\\bar{t}$ ",normed=True,bins=100)
# lg = plt.legend()
# lg.draw_frame(False)
# plt.gca().axes.xaxis.set_ticklabels([])
# #hide the 0
# plt.gcf().canvas.draw()
# A=plt.gca().axes.yaxis.get_ticklabels()
# B=[y.get_text() for y in A]
# B[0]=u""
# plt.gca().axes.yaxis.set_ticklabels(B)
# ax2=plt.subplot(gs2[0, 0]) # now we plot
# ax2.tick_params(axis='x', labelsize=25)
# ax2.tick_params(axis='y', labelsize=10)
# # ax2.tick_params(axis='y', labelsize=25)
# plt.hist(d0Bs[abs(d0Bs)>0],color="black",range=[-8.0,8.0],histtype="step",lw=2,label=" $t\\bar{t}$ ",normed=True,bins=100)
# plt.ylim([0,0.01])
# plt.xlabel(' $|d_{0}|$ of tracks coming from b vertex [mm]',fontsize=30)
# plt.show()



plt.show()



