import glob
import ROOT
import os
import math


execfile("../plot_helpers/basic_plotting.py")
execfile("signalMap.py")

def firstLep(event, i_ignore):
    #print i_ignore 
    l0 = -1
    pt0 = -1

    for il in xrange(len(event.truthLepton_pt)):
        if event.truthLepton_pt[il] > pt0 and (abs(event.truthLepton_pdgId[il]) == 11 or abs(event.truthLepton_pdgId[il]) == 13) and il!= i_ignore:
            pt0 = event.truthLepton_pt[il]
            l0 = il

    return l0


reco_weight_t = "lumi_weight*mcEventWeight*muSF*mu_sel_sf*mu_trig_sf*disp_sf*el_reco_sf*pileupWeight"

d = "/eos/atlas/atlascerngroupdisk/phys-susy/Highd0Lep_ANA-SUSY-2018-14/stat_ext_new/signal/"

append = "slep"
if 'slep' in append:
    lifetimes = ["0p01", "0p1", "1", "10"]
    masses = ["50","100","200",  "300", "400", "500","600","700","800","900","1000"]
else:
    lifetimes = ["0p01", "0p1", "1"]
    masses = ["50","100","200",  "300", "400", "500"]


SRs = ["all","ee","em","mm"]


# Set up hists for acceptance
hists = {}
for sr in SRs:
    hists[sr] = {}

    if 'slep' in append:
        hists[sr]["accepted"] =     ROOT.TH2F(sr+"_accepted", sr+"_accepted",   21, 0, 1050, 5, -2., 3.)
        hists[sr]["selected"] =    ROOT.TH2F(sr+"_selected", sr+"_selected", 21, 0, 1050, 5, -2., 3.)
        hists[sr]["total"] =     ROOT.TH2F(sr+"_total", sr+"_total",   21, 0, 1050, 5, -2., 3.)
    else:
        hists[sr]["accepted"] =     ROOT.TH2F(sr+"_pass", sr+"_pass",   11, 0, 550, 4, -2., 2.)
        hists[sr]["selected"] =    ROOT.TH2F(sr+"_total", sr+"_total", 11, 0, 550, 4, -2., 2.)
        hists[sr]["total"] =     ROOT.TH2F(sr+"_truth", sr+"_truth", 11, 0, 550, 4, -2., 2.)


for lt in lifetimes:
    for m in masses:

        print "Working with mass", m, "and lifetime", lt

        mass = int(m)
        if lt == "0p01": lifetime = -2
        elif lt == "0p1": lifetime = -1
        elif lt == "1": lifetime = 0
        elif lt == "10": lifetime = 1
        
        slep_str = "SlepSlep_"
        if "stau" in append: slep_str = "StauStau_"
        sig_str = slep_str + m + "_0_" + lt + "ns"
        print sig_str

         
        try:
            dsid = signals_inv[sig_str]
        except:
            print "no dsid found for ", sig_str
            continue

        fstr = "%s%s*.root"%(d,dsid)
        t = getTree(fstr)  
        if t==-1: continue
        print "Found events:", t.GetEntries()

        maxEntries = t.GetEntries()
        
        filename = glob.glob(fstr)[0]
        xs = getXS(dsid)
        print "cross section: ", xs
        
            
        for event in t:
            
            reco_weight = event.mcEventWeight*event.muSF*event.mu_sel_sf*event.mu_trig_sf*event.disp_sf*event.el_reco_sf*event.pileupWeight*event.lumi_weight/(139000.*xs)
            

            for s in SRs :
                hists[s]["total"].Fill(mass,lifetime)

            ##### efficiency numerator 
            sr = ""
            if event.SRmm: sr = "mm"
            if event.SRee: sr = "ee" 
            if event.SRem: sr = "em"
            if sr != "":
                hists[sr]["selected"].Fill(mass,lifetime,reco_weight)
                hists["all"]["selected"].Fill(mass,lifetime,reco_weight)
                
                m0 = getTruthMatch(event.lepton_truthMatchedBarcode[0])
                m1 = getTruthMatch(event.lepton_truthMatchedBarcode[1])
                if m0 != -1 and m1 != -1:
                    hists[sr]["pt_d0_selected_lead"].Fill(event.truthLepton_pt[m0], abs(event.truthLepton_d0[m0]))
                    hists[sr]["pt_d0_selected_sublead"].Fill(event.truthLepton_pt[m1], abs(event.truthLepton_d0[m1]))
                
             
                

            ## evens must have 2 leptons
            if len(event.truthLepton_pt) < 2 : continue
            
            ## define signal regions based on the two highest pt electrons/muons
            il0 = firstLep(event, -1) 
            il1 = firstLep(event, il0) 
            
            if il0 == -1 or il1 == -1: continue
           
            ## in practice, there is only one different flavor region, SR-em (which will be collapsed and is only separated to make the logic in this script easier) 
            sr = ""
            if abs(event.truthLepton_pdgId[il0]) == 11 and abs(event.truthLepton_pdgId[il1]) == 11: sr = "ee"
            if abs(event.truthLepton_pdgId[il0]) == 13 and abs(event.truthLepton_pdgId[il1]) == 13: sr = "mm"
            if abs(event.truthLepton_pdgId[il0]) == 11 and abs(event.truthLepton_pdgId[il1]) == 13: sr = "em"
            if abs(event.truthLepton_pdgId[il0]) == 13 and abs(event.truthLepton_pdgId[il1]) == 11: sr = "me"
            
            
            ## collapse SRem and SRme for histogram filling
            srf = sr
            if sr == "me": srf = "em"
            
            
            ##### now calculate acceptance
            
            # Different eta cuts for electrons and muons -- based on atlas recommendations
            if (sr == "ee" or sr == "em") and abs(event.truthLepton_eta[il0]) > 2.47: continue 
            if (sr == "ee" or sr == "me") and abs(event.truthLepton_eta[il1]) > 2.47: continue 
            if (sr == "mm" or sr == "me") and abs(event.truthLepton_eta[il0]) > 2.5: continue 
            if (sr == "mm" or sr == "em") and abs(event.truthLepton_eta[il1]) > 2.5: continue 
            

            
            ## signal pt and d0 cuts 
            if event.truthLepton_pt[il0] < 65 or event.truthLepton_pt[il1] < 65: continue

            if abs(event.truthLepton_d0[il0]) < 3 or abs(event.truthLepton_d0[il1]) < 3: continue

           

            ## triggers are applied depending on the kinematics of the event
            ##   the pt cuts are higher than those in the trigger to match the cuts used in 
            ##   the filters that determine which events get reconstructed with special tracking
            ##                  - HLT_g140_loose if there is an electron with pt > 160 GeV 
            ##                  - HLT_2g50_loose if there are 2 electrons each with pt > 60 GeV 
            ##                  - HLT_m60_0eta105_msonly if there is a muon with pt > 60 GeV and |eta| < 1.07
            ## Events must fall into one of these regions to be accepted, and then the appropriate trigger must pass for the event to be selected 
            tr = False
            if (sr == "ee" or sr == "em") and event.truthLepton_pt[il0] > 160: tr = True
            if (sr == "ee" or sr == "me") and event.truthLepton_pt[il1] > 160: tr = True
            if sr == "ee" and event.truthLepton_pt[il0] > 60 and event.truthLepton_pt[il1] > 60: tr = True
            if (sr == "mm" or sr == "me") and (event.truthLepton_pt[il0] > 60 and abs(event.truthLepton_eta[il0]) < 1.07): tr = True
            if (sr == "mm" or sr == "em") and (event.truthLepton_pt[il1] > 60 and abs(event.truthLepton_eta[il1]) < 1.07): tr = True
            if not tr: continue 
            


            ## delta R cut between two leptons to minimize material interactions
            lep0 = ROOT.TVector3()
            lep0.SetPtEtaPhi(event.truthLepton_pt[il0],event.truthLepton_eta[il0], event.truthLepton_phi[il0])
            lep1 = ROOT.TVector3()
            lep1.SetPtEtaPhi(event.truthLepton_pt[il1], event.truthLepton_eta[il1], event.truthLepton_phi[il1])
            
            if abs(lep0.DeltaR(lep1)) < 0.2: continue

            
            ## acceptance numerator/reco denominator 
            hists[srf]["accepted"].Fill(mass,lifetime, event.mcEventWeight)
            hists["all"]["accepted"].Fill(mass,lifetime, event.mcEventWeight)
            
            


