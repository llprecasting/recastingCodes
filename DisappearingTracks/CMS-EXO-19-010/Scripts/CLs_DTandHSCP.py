import numpy as np
import random
import math

import os

import math
import sys

import re
from collections import OrderedDict

#from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import interp1d
#from scipy import interpolate

#from math import log10

try:
    import scipy.stats
except ImportError:
    print("Cannot import scipy")
    raise


from statistics import mean



def Proportions(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
       

    ### NumToyExperiments is a maximum. We'll keep going until we have an accurate enough answer.
    
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    min_proportional_error=0.025
    min_absolute_error=0.005
    properrB=1.0
    properrSplusB=1.0
    pts=0
    batch_size=5000
    pbs = []
    
    pSplusBs = []
    totalB=0
    totalBplusS=0
    ErrorTooLarge=True
    ## proportional error on background, absolute error on signal
    while ErrorTooLarge and pts < NumToyExperiments:
        
        pts = pts + batch_size
        ExpectedBGs = scipy.stats.norm.rvs(loc=ExpectedBG, scale=BGError, size=batch_size)

        # Ignore values in the tail of the Gaussian extending to negative numbers
        ExpectedBGs = [value for value in ExpectedBGs if value > 0]
    
        # For each toy experiment, get the actual number of background events by
        # taking one value from a Poisson distribution created using the expected
        # number of events.
        ToyBGs = scipy.stats.poisson.rvs(ExpectedBGs)
        ToyBGs = list(map(float, ToyBGs))

        # The probability for the background alone to fluctutate as LOW as
        # observed = the fraction of the toy experiments with backgrounds as low as
        # observed = p_b.
        # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
        #print(ToyBGs)
        #print(NumObserved)
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.percentileofscore.html
        # 'weak': This kind corresponds to the definition of a cumulative distribution function. A percentileofscore of 80% means that 80% of values are less than or equal to the provided score.
        #p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01
        batch_p_b=scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01
        totalB=totalB+len(ToyBGs)
        #print("pb: %.3f" %p_b)
        # Toy MC for background+signal
        ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
        ExpectedBGandS = [x for x in ExpectedBGandS if x > 0]
        
        if len(ExpectedBGandS)==0:
            continue
            #return 0., p_b
        ToyBplusS = scipy.stats.poisson.rvs(ExpectedBGandS)
        ToyBplusS = list(map(float, ToyBplusS))
        totalBplusS=totalBplusS+len(ToyBplusS)
        
        # Calculate the fraction of these that are >= the number observed,
        # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
        batch_p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01
        pSplusBs.append(batch_p_SplusB*len(ToyBplusS)/batch_size)
        pbs.append(batch_p_b*len(ToyBGs)/batch_size)
        p_b=np.mean(pbs)
        p_SplusB=np.mean(pSplusBs)

        ErrorTooLarge=False
        ### Proportional error on background, because if p_b is small need to be quite sure of it ... 
        if p_b > 0.1:
            ## Accept absolute error
            properrB=math.sqrt(p_b*(1-p_b)/(totalB))
            if properrB > min_absolute_error:
                ErrorTooLarge=True
        elif p_b > 0.0:
            ## Require proportional error
            properrB=math.sqrt((1-p_b)/(p_b*totalB))
            if properrB > min_proportional_error:
                ErrorTooLarge=True
            #properrB=math.sqrt(p_b*(1-p_b)/(totalB))
        else:
            properrB=0.0
            ## Not sure what to do about the error here ... keep going?
            ErrorTooLarge=True
            
        ## Not proportional error anymore, just absolute error on the probability. Actually that is all we care about
        # if p_SplusB > 0.0:
        #     #properrSplusB=math.sqrt((1-p_SplusB)/(p_SplusB*totalBplusS))
        #     properrSplusB=math.sqrt(p_SplusB*(1-p_SplusB)/(totalBplusS))
        # else:
        #     properrSplusB=0.0

        # Check for uncertainty on signal, 
        # but also check we are not wasting our time on very low or very high probability points
        
        if p_SplusB > 0.1:
            ## Accept absolute error
            properrSplusB=math.sqrt(p_SplusB*(1-p_SplusB)/(totalBplusS))
            if properrSplusB > min_absolute_error:
                ErrorTooLarge=True
            ## Don't waste time on points with small signal
            if p_b > 0.1:
                if (p_SplusB-3*properrSplusB)/p_b > 0.5:
                    ErrorTooLarge=False 
        elif p_SplusB > 0.0:
            ## Require proportional error
            properrSplusB=math.sqrt((1-p_SplusB)/(p_SplusB*totalBplusS))
            if properrSplusB > min_proportional_error:
                ErrorTooLarge=True
            #properrB=math.sqrt(p_b*(1-p_b)/(totalB))
            ## Don't waste time on points with very large signal!
            if p_b > 0.1:
                if p_SplusB*(1+3*properrSplusB)/p_b < 0.01:
                    ErrorTooLarge=False 
        else:
            properrSplusB=0.0
            ## Presumably this is ok that we have a very small probability of signal
            #ErrorTooLarge=True
        
        
        
        
        #err=max(properrB,properrSplusB)
    #print("%d %.4f (%.4f) %.4f (%.4f)" %(pts,p_b,properrB,p_SplusB,properrSplusB))
    return p_SplusB, p_b


def CLs(Observed, ExpectedBGs, BGErrors, SigHypotheses, NumToyExperiments):
    p_sb=1.0
    p_b=1.0

    for x,y,z,u in zip(Observed,ExpectedBGs,BGErrors,SigHypotheses):
        #print(x,y,z,u)
        try:
            tp_sb,tp_b = Proportions(x,y,z,u,NumToyExperiments)
            if tp_sb > tp_b:
                tCLs = 0.0
            else:
                tCLs = 1.0 - tp_sb/tp_b
            #print("s+b: %.3f b: %.3f, region CLs: %.3f" %(tp_sb,tp_b,tCLs))
        except:
            print("failed")
            return False

        p_sb=p_sb*tp_sb
        p_b=p_b*tp_b

    if p_sb>p_b:
        #print(" CLs: %.3f" %(0.0))
        return 0.
    else:
        #print(" CLs: %.3f" %(1.-(p_sb / p_b)))
        return 1.-(p_sb / p_b) # 1 - CLs





def GetLimits(Observed,backgrounds,BGerrors,signals):
    ## First get the raw CLs
    temptoys=500000
    mycls=CLs(Observed,backgrounds,BGerrors,signals,temptoys)

    return mycls

def GetSigLimit(Observed,backgrounds,BGerrors,refSigs):

    #refxs=XSdata_chch.getxs(pointmass)+XSdata_chnu.getxs(pointmass)
    temptoys=100000
    if max(refSigs) ==0.0:
        raise
    def GetSig95(x):  
        nsigs=[x * y for y in refSigs]
        if max(nsigs) < 0.001:
            return -0.95
        return CLs(Observed,backgrounds,BGerrors,nsigs,temptoys)-0.95

    mlow = 1
    mhig = 1

    while GetSig95(mlow)>0.0:
        mlow=mlow*0.1
 
    while GetSig95(mhig)<0.0:
        mhig=mhig*10.0

    #print("Testing between %.3f and %.3f " %(mlow,mhig))
    try:
        m95 = scipy.optimize.brentq(GetSig95,mlow,mhig,xtol=mlow/100.)
    except:
        m95=-1

    return m95

#### DATA FOR THE CMS DISAPPEARING TRACK ANALYSIS


lumis=OrderedDict()




lumis['SR3_2015']=2.7
#lumis['SR3_2016']=35.7

lumis['SR3_2016A']=8.3
lumis['SR3_2016B']=27.4


lumis['SR1_2017']=42.0
lumis['SR1_2018A']=21.0
lumis['SR1_2018B']=39.0

lumis['SR2_2017']=42.0
lumis['SR2_2018A']=21.0
lumis['SR2_2018B']=39.0

lumis['SR3_2017']=42.0
lumis['SR3_2018A']=21.0
lumis['SR3_2018B']=39.0
              
backgrounds=OrderedDict()

backgrounds['SR3_2015']=[0.1,0.1]
#backgrounds['SR3_2016']=[6.4,1.4]
backgrounds['SR3_2016A']=[2.4,0.6]
backgrounds['SR3_2016B']=[4.0,1.1]


backgrounds['SR1_2017']=[12.2,4.8]
backgrounds['SR1_2018A']=[7.3,3.5]
backgrounds['SR1_2018B']=[10.3,5.4]

backgrounds['SR2_2017']=[2.1,0.7]
backgrounds['SR2_2018A']=[0.6,0.3]
backgrounds['SR2_2018B']=[1.0,0.3]

backgrounds['SR3_2017']=[6.7,1.3]
backgrounds['SR3_2018A']=[1.8,0.6]
backgrounds['SR3_2018B']=[5.7,1.3]


##### DATA FOR THE HSCP ATLAS ANALYSIS

cand2b=[1.5,0.06,0.007,0.0017]
cand2bu=[0.3,0.01,0.002,0.0009]
cand2o=[0,0,0,0]

cand1b=[240.0,17.0,2.2,0.48]
cand1bu=[20.0,2.0,0.2,0.07]
cand1o=[227,16,1,0]

AllRegions=OrderedDict()

AllRegions['DT_CMS']=[['SR1_2017','SR1_2018A','SR1_2018B'],['SR2_2017','SR2_2018A','SR2_2018B'],['SR3_2015','SR3_2016A','SR3_2016B','SR3_2017','SR3_2018A','SR3_2018B']]

AllRegions['HSCP_ATLAS']=[]

#regionlists=[['SR1_2017','SR1_2018A','SR1_2018B'],['SR2_2017','SR2_2018A','SR2_2018B'],['SR3_2015','SR3_2016A','SR3_2016B','SR3_2017','SR3_2018A','SR3_2018B'],['CAND1_0'],['CAND1_1'],['CAND1_2'],['CAND1_3'],['CAND2_0'],['CAND2_1'],['CAND2_2'],['CAND2_3']]


#observed={'SR3_2015':1, 'SR3_2016A':2, 'SR3_2016B':4, 'SR1_2017':17,'SR1_2018A':5,'SR1_2018B':11,'SR2_2017':4,'SR2_2018A':0,'SR2_2018B':2,'SR3_2017':6,'SR3_2018A':2,'SR3_2018B':1, 'CAND1_0': 227, 'CAND1_1':16, 'CAND1_2':1,'CAND_1_3':0,'CAND2_0':0,'CAND2_1':0,'CAND2_2':0,'CAND2_3':0}
observed={'SR3_2015':1, 'SR3_2016A':2, 'SR3_2016B':4, 'SR1_2017':17,'SR1_2018A':5,'SR1_2018B':11,'SR2_2017':4,'SR2_2018A':0,'SR2_2018B':2,'SR3_2017':6,'SR3_2018A':2,'SR3_2018B':1}

#regionlists=[['SR1_2017','SR1_2018A','SR1_2018B'],['SR2_2017','SR2_2018A','SR2_2018B'],['SR3_2015','SR3_2016A','SR3_2016B','SR3_2017','SR3_2018A','SR3_2018B'],['CAND1_0'],['CAND1_1'],['CAND1_2'],['CAND1_3'],['CAND2_0'],['CAND2_1'],['CAND2_2'],['CAND2_3']]
#regionlists=[['SR1_2017','SR1_2018A','SR1_2018B'],['SR2_2017','SR2_2018A','SR2_2018B'],['SR3_2015','SR3_2016A','SR3_2016B','SR3_2017','SR3_2018A','SR3_2018B']]

for k in [0,1,2,3]:
    c1s='CAND1_'+str(k)
    c2s='CAND2_'+str(k)
    lumis[c1s]=36.1
    lumis[c2s]=36.1
    backgrounds[c1s]=[cand1b[k],cand1bu[k]]
    backgrounds[c2s]=[cand2b[k],cand2bu[k]]
    observed[c1s]=cand1o[k]
    observed[c2s]=cand2o[k]
    AllRegions['HSCP_ATLAS'].append([c1s])
    AllRegions['HSCP_ATLAS'].append([c2s])
    #regionlists.append([c1s])
    #regionlists.append([c2s])



def ReadMyFile(fname):
    ### Need to potentially read two sets of data
    ### Only interested in exporting efficiencies
    EffBlocks=OrderedDict()
    EffBlock=OrderedDict()
    
    try:
        IF=open(fname,'r')
    except:
        raise
    ReadingEff= False
    ReadingProcess = False
    ANAME=""
    for line in IF:
        if line.startswith('==='): ## reached end of useful bit
            break
        sline=line.strip().split('#')[0]
        sline=line.replace('\t',' ')
        #sline=line.replace('  ',' ')
        parts=sline.split()
        if len(parts) < 2:
            continue
        if sline.upper()[0:6]=='BLOCK ':
        ## append old block
            blockname=parts[1].upper()
            if blockname == "EFFICIENCIES":
                ReadingEff=True
                ReadingProcess=False
            elif blockname == "PROCESS":
                
                #print("Found new process")
                if ReadingProcess or ReadingEff:  ### already have at least one, so add to stack
                    if ANAME=="":
                        # try to guess the analysis name based on the signal regions
                        if 'SR1_2017' in EffBlock:
                            ANAME='DT_CMS'
                        elif 'CAND1_0' in EffBlock:
                            ANAME='HSCP_ATLAS'
                        else:
                            ANAME='UNKNOWN'
                    if ANAME=='DT_CMS' and 'SR3_2016' in EffBlock:
                        EffBlock['SR3_2016A']=EffBlock['SR3_2016']
                        EffBlock['SR3_2016B']=EffBlock['SR3_2016']
                    EffBlocks[ANAME]=EffBlock
                    if len(parts) > 2:
                        ANAME=parts[2]
                    else:
                        ANAME=""     
                    #EffBlocks.append(EffBlock)
                    EffBlock=OrderedDict()
                #print(EffBlock)
                ReadingEff=False
                ReadingProcess=True
                
            else:
                ReadingEff=False
                ReadingProcess=False
        else:
            if ReadingEff or ReadingProcess:
                EffBlock[parts[0]]=float(parts[1]) ## don't care about errors
    #EffBlocks.append(EffBlock)
    if ANAME=="":
        ## try to guess the analysis name based on the signal regions
        if 'SR1_2017' in EffBlock:
            ANAME='DT_CMS'
        elif 'CAND1_0' in EffBlock:
            ANAME='HSCP_ATLAS'
        else:
            ANAME='UNKNOWN'
    if ANAME=='DT_CMS' and 'SR3_2016' in EffBlock:
        EffBlock['SR3_2016A']=EffBlock['SR3_2016']
        EffBlock['SR3_2016B']=EffBlock['SR3_2016']
    EffBlocks[ANAME]=EffBlock
    #print(EffBlocks)
    IF.close()  
    return EffBlocks

def GetLimitsFromFile(fname,inXS=False):
    #print('Starting limits with '+fname)
    AllData=OrderedDict()
    try:
        AllData=ReadMyFile(fname)
    except:
        print('Failed to read file')
        raise
    #print(AllData)
    for analysis in AllData:
        print('Getting limits for '+analysis)
        AnalysisData=OrderedDict()
        myeffs=AllData[analysis]
        if not inXS:
            XS=myeffs['XS']
        else:
            XS=inXS

        allsigs=[]
        allbacks=[]
        alluncerts=[]
        allobserveds=[]
        
        obsCLsvals=[]
        expCLsvals=[]
        lim95s=[]
        bestgroup=None
        regionlists=AllRegions[analysis]
        #print(regionlists)
        for k,tocombine in enumerate(regionlists):
            print('Region %d is %s' %(k,tocombine))
            sigs=[]
            backs=[]
            uncerts=[]
            observeds=[]
            for region in tocombine:
                backs.append(backgrounds[region][0])
                uncerts.append(backgrounds[region][1])
                observeds.append(observed[region])
                lumi=lumis[region]
                myn=myeffs[region]*lumi
                sigs.append(myn*XS)
            
            allbacks.append(backs)
            alluncerts.append(uncerts)
            allobserveds.append(observeds)
            allsigs.append(sigs)
            ### now have the signals, backgrounds etc for the analysis
            #print(sigs)
            expCLs=GetLimits(backs,backs,uncerts,sigs)
            
            obsCLs=GetLimits(observeds,backs,uncerts,sigs)
            #print('%.3e %.3e'%(expCLs,obsCLs))
            obsCLsvals.append(obsCLs)
            expCLsvals.append(expCLs)
        
        bestexp=-1.0
        bestobs=-1.0
        bestreg=-1
        for k,exp in enumerate(expCLsvals):
            if exp > bestexp:
                bestreg=k
                bestexp=exp
                bestobs=obsCLsvals[k]
        AnalysisData['Best Region']=bestreg
        AnalysisData['Best CLs']=bestobs
        AnalysisData['All Expected']=expCLsvals
        AnalysisData['All Observed']=obsCLsvals
        #print('here')
        if bestreg > -1:
            k=bestreg
            backs=allbacks[k]
            uncerts=alluncerts[k]
            observeds=allobserveds[k]
            sigs=allsigs[k]
            ## only bother to find the upper limit on
            try:
                sig95=XS*GetSigLimit(observeds,backs,uncerts,sigs)
            except:
                sig95=-1
            AnalysisData['s95']=sig95
        AllData[analysis]=AnalysisData

    return AllData



            
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Need to specify an efficiency file name!')
        print('Usage: '+sys.argv[0]+' <efficiency filename> [optional: cross-section]')
        raise SystemExit
    
    fname=sys.argv[1]
    if len(sys.argv) > 2:
        try:
            inXS=float(sys.argv[2])
        except:
            print('Input cross-section must be a number!')
            print('Usage: '+sys.argv[0]+' <efficiency filename> [optional: cross-section]')
            raise SystemExit
        try:
            alldata=GetLimitsFromFile(fname,inXS)
        except:
            print('Failed to get limits')
            raise SystemExit

    else:
        try:
            alldata=GetLimitsFromFile(fname)
        except:
            print('Failed to get limits')
            raise SystemExit
 
    for analysis in alldata:
        print('Analysis: %s, Best region: %d, Best CLs limit: %.3f, Upper cross-section limit: %.3e fb' %(analysis,alldata[analysis]['Best Region'],alldata[analysis]['Best CLs'],alldata[analysis]['s95']))
