#################################################
# Plotting script for displaced vertex recast
# Written by Giovanna Cottin (gfcottin@gmail.com)
#################################################
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
##########################
#Data for ATLAS Benchmarks
##########################
#Eff Vs rDV
ATLAS_effVsrDV_file = "ATLAS_data/rDV_eff.dat"
ATLAS_effVsrDV_data = {0:[],1:[]}
for line in open(ATLAS_effVsrDV_file,"r"):
  for i,C in enumerate(line.split()): ATLAS_effVsrDV_data[i].append(float(C))
ATLAS_rDV    =  np.array(ATLAS_effVsrDV_data[0])
ATLAS_effrDV =  np.array(ATLAS_effVsrDV_data[1])
# ATLAS Split SUSY mchi = 100 GeV ; m_gluino = 1.4 TeV
# ATLAS_gluino1_file  = "ATLAS_data/limit_plot1_1400.dat"
ATLAS_gluino1_file  = "ATLAS_data/limits_gluino1_NEW.dat"
ATLAS_gluino1_data  = {0:[],1:[]}
for line in open(ATLAS_gluino1_file,"r"):
  for i,C in enumerate(line.split()): ATLAS_gluino1_data[i].append(float(C))
gluino1_tau      = np.array(ATLAS_gluino1_data[0])
gluino1_crossSec = np.array(ATLAS_gluino1_data[1])*1000.#in [fb]
# ATLAS Split SUSY mchi = 100 GeV ; m_gluino = 2 TeV
# ATLAS_gluino2_file  = "ATLAS_data/limit_plot1_2000.dat"
ATLAS_gluino2_file  = "ATLAS_data/limits_gluino2_NEW.dat"

ATLAS_gluino2_data  = {0:[],1:[]}
for line in open(ATLAS_gluino2_file,"r"):
  for i,C in enumerate(line.split()): ATLAS_gluino2_data[i].append(float(C))
gluino2_tau      = np.array(ATLAS_gluino2_data[0])
gluino2_crossSec = np.array(ATLAS_gluino2_data[1])*1000.#in [fb]
##################
#My eff data files
##################
effVsLifetime  = "truthDV_data/effVsTau_1400_RhadID_FIXED.dat"

effVsLifetime_data  = {0:[],1:[],2:[]}
for line in open(effVsLifetime,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_data[i].append(float(C))
tau      = np.array(effVsLifetime_data[0])/299.792 # in ns
nPass    = np.array(effVsLifetime_data[1])
eff      = np.array(effVsLifetime_data[2])
deltaEff = 1.0/np.sqrt(nPass)
effVsLifetime2  = "truthDV_data/effVsTau_2000_RhadID_FIXED.dat"
effVsLifetime_data2  = {0:[],1:[],2:[]}
for line in open(effVsLifetime2,"r"):
  for i,C in enumerate(line.split()): effVsLifetime_data2[i].append(float(C))
tau2      = np.array(effVsLifetime_data2[0])/299.792## in ns
nPass2    = np.array(effVsLifetime_data2[1])
eff2      = np.array(effVsLifetime_data2[2])
deltaEff2 = 1.0/np.sqrt(nPass2)

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
          'legend.fontsize': 15
          # 'legend.fontsize': 20
          }
plt.rcParams.update(params) 
###########################
lp=0.27
rp=0.75
tp=0.94
bp=0.14
#modify bottom y left
#for n subplots is /(n+(n-1)wp  
# ####################################
plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig1=plt.figure(1,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3) 
gs2.update(top=0.3,bottom=bp)  
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)
plt.ylim([0.0, 1.0])
plt.xscale('log')
ax1.text(0.0015,0.9 , "split SUSY simplified model", fontsize=15)
ax1.text(0.0015,0.8 , "$\\tilde{g} \\rightarrow q \\bar{q} \\tilde{\\chi}^{0}_{1}$", fontsize=15)
ax1.text(0.0015,0.6 , "$m_{\\tilde{g}} = 1400 $ GeV", fontsize=15)
ax1.text(0.0015,0.5 , "$m_{\\tilde{\\chi}^{0}_{1}} = 100$ GeV", fontsize=15)
plt.plot(tau, eff, color='red',linestyle='--',lw=2,label=" Recasted")
plt.fill_between(tau,eff-deltaEff*eff,eff+deltaEff*eff,color='red',alpha=0.3)
lg = plt.legend(loc=1)
lg.draw_frame(False)
plt.ylabel('event efficiency',fontsize=15)
plt.xlabel('$\\tau$ [ns]',fontsize=15)
plt.savefig('Plots/effVsTau_1400_FIXED.pdf', bbox_inches='tight')


plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig11=plt.figure(11,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3) 
gs2.update(top=0.3,bottom=bp)  
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)
plt.ylim([0.0, 1.0])
plt.xscale('log')
ax1.text(0.0015,0.9 , "split SUSY simplified model", fontsize=15)
ax1.text(0.0015,0.8 , "$\\tilde{g} \\rightarrow q \\bar{q} \\tilde{\\chi}^{0}_{1}$", fontsize=15)
ax1.text(0.0015,0.6 , "$m_{\\tilde{g}} = 2000 $ GeV", fontsize=15)
ax1.text(0.0015,0.5 , "$m_{\\tilde{\\chi}^{0}_{1}} = 100$ GeV", fontsize=15)
plt.plot(tau2, eff2, color='red',linestyle='--',lw=2,label=" Recasted")
plt.fill_between(tau2,eff2-deltaEff2*eff2,eff2+deltaEff2*eff2,color='red',alpha=0.3)
lg = plt.legend(loc=1)
lg.draw_frame(False)
plt.ylabel('event efficiency',fontsize=15)
plt.xlabel('$\\tau$ [ns]',fontsize=15)
plt.savefig('Plots/effVsTau_2000_FIXED.pdf', bbox_inches='tight')
#####################################################################################
### Approx upper limit on sigma (3-0.02)/32.7
### 0.09113149847094801
limit1=0.0911314/eff
limit2=0.0911314/eff2
# eff[0] corresponds to largest tau, gluino_crossSec starts with lowest tau. I need a consistent difference
diff1 = gluino1_crossSec/limit1[::-1] 
diff2 = gluino2_crossSec/limit2[::-1] 


plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig3=plt.figure(3,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3) 
gs2.update(top=0.3,bottom=bp)  
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.text(0.005,50 , "split SUSY simplified model", fontsize=15)
ax1.text(0.005,30 , "$\\tilde{g} \\rightarrow q \\bar{q} \\tilde{\\chi}^{0}_{1}$", fontsize=15)
ax1.text(0.005,15 , "$m_{\\tilde{g}} = 1400 $ GeV", fontsize=15)
ax1.text(0.005,10 , "$m_{\\tilde{\\chi}^{0}_{1}} = 100$ GeV", fontsize=15)
plt.xscale('log')
plt.yscale('log')
plt.plot(gluino1_tau,gluino1_crossSec,color="black",lw=2,label="ATLAS")
plt.plot(tau, limit1, color='red',linestyle='--',lw=2, label=" Recasted")
plt.fill_between(tau,limit1-deltaEff*limit1,limit1+deltaEff*limit1,color='red',alpha=0.3)
lg=plt.legend(scatterpoints=1)
lg.draw_frame(False)
# plt.xlabel('$\\tau$ [ns]',fontsize=20)
plt.ylabel('upper limit on cross section [fb]',fontsize=20)
plt.gca().axes.xaxis.set_ticklabels([])
## hide the 0
plt.gcf().canvas.draw()
A=plt.gca().axes.yaxis.get_ticklabels()
B=[y.get_text() for y in A]
B[0]=u""
plt.gca().axes.yaxis.set_ticklabels(B)
ax2=plt.subplot(gs2[0, 0]) # now we plot
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=10)
plt.plot(gluino1_tau,diff1,color='red',lw=2)
plt.axhline(y=1,color='k',linestyle='--',lw=2,alpha=0.3)
# plt.ylim(-1.0,1.0)
# plt.ylim(-2.0,2.0)
plt.ylim(0.5,2.0)
plt.xscale('log')
plt.ylabel('ratio',fontsize=15)
plt.xlabel('$\\tau$ [ns]',fontsize=20)
plt.savefig('Plots/limits_gluino1_FIXED.pdf', bbox_inches='tight')

plotaspect=(1-lp-(1-rp))/(1-bp-(1-tp)) #P=(1-l-r)/(1-t-b)/(n+(n-1)*w)*(n+(n-1)*h)
fig4=plt.figure(4,figsize=(1.5*(3.37*2), 1.5*(3.37*2*plotaspect)))
plt.subplots_adjust(left=lp,right=rp,top=tp,bottom=bp)
gs1 = gridspec.GridSpec(1, 1) 
gs2 = gridspec.GridSpec(1, 1)
gs1.update(top=tp, bottom=0.3) 
gs2.update(top=0.3,bottom=bp)  
ax1=plt.subplot(gs1[0, 0])
ax1.tick_params(axis='y', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.text(0.009,250 , "split SUSY simplified model", fontsize=15)
ax1.text(0.009,150 , "$\\tilde{g} \\rightarrow q \\bar{q} \\tilde{\\chi}^{0}_{1}$", fontsize=15)
ax1.text(0.009,60 , "$m_{\\tilde{g}} = 2000 $ GeV", fontsize=15)
ax1.text(0.009,35 , "$m_{\\tilde{\\chi}^{0}_{1}} = 100$ GeV", fontsize=15)
plt.xscale('log')
plt.yscale('log')
plt.plot(gluino2_tau,gluino2_crossSec,color="black",lw=2,label="ATLAS")
plt.plot(tau2, limit2, color='red',linestyle='--',lw=2, label=" Recasted")
plt.fill_between(tau2,limit2-deltaEff2*limit2,limit2+deltaEff2*limit2,color='red',alpha=0.3)
lg=plt.legend(scatterpoints=1)
lg.draw_frame(False)
#plt.xlabel('$\\tau$ [ns]',fontsize=20)
plt.ylabel('upper limit on cross section [fb]',fontsize=20)
plt.gca().axes.xaxis.set_ticklabels([])
#hide the 0
plt.gcf().canvas.draw()
A=plt.gca().axes.yaxis.get_ticklabels()
B=[y.get_text() for y in A]
B[0]=u""
plt.gca().axes.yaxis.set_ticklabels(B)
ax2=plt.subplot(gs2[0, 0]) # now we plot
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=10)
plt.plot(gluino2_tau,diff2,color='red',lw=2)
plt.axhline(y=1,color='k',linestyle='--',lw=2,alpha=0.3)
# plt.ylim(-1.0,1.0)
plt.ylim(0.5,2.0)
plt.xscale('log')
plt.ylabel('ratio',fontsize=15)
plt.xlabel('$\\tau$ [ns]',fontsize=20)
plt.savefig('Plots/limits_gluino2_FIXED.pdf', bbox_inches='tight')
###################################
plt.show()
