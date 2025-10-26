#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import glob
import re
import pickle

def runmatcher(patt, filepath):
    return float(re.search(patt,filepath)[1])

class normalizer:
    def __init__(self, data_dir ='.'):
        self.pp2chgchg_path = Path(data_dir) / 'xsec_13TeV_pp2chgchg_fb.csv'
        self.pp2chgn1_path  = Path(data_dir) / 'xsec_13TeV_pp2chgn1_fb.csv'
        self.pp2chgchg_ref = pd.read_csv(self.pp2chgchg_path,header=0, sep=',')
        self.pp2chgn1_ref  = pd.read_csv(self.pp2chgn1_path ,header=0, sep=',')

    def xsec_13TeV_fb(self, mass):
        return {'pp2chgchg': np.interp(mass,self.pp2chgchg_ref['mass'].to_numpy(),self.pp2chgchg_ref['xsec'].to_numpy()),'pp2chgn1': np.interp(mass,self.pp2chgn1_ref['mass'].to_numpy(),self.pp2chgn1_ref['xsec'].to_numpy())}

    def ratios(self,mass):
        ref_xsecs = self.xsec_13TeV_fb(mass)
        sum = ref_xsecs['pp2chgchg'] + ref_xsecs['pp2chgn1']
        return {key:value/sum for key,value in ref_xsecs.items()}

    def nc(self, mass, val_chgchg, val_chgn1):
        chn_ratios = self.ratios(mass)
        return val_chgchg*ratios['pp2chgchg'] + val_chgn1*ratios['pp2chgn1']

    def normalize_evts(self, mass, SRevts, totalevts):
        ref_xsecs = self.xsec_13TeV_fb(mass)
        return {key:value*137*SRevts[key]/totalevts[key] for key,value in ref_xsecs.items()}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate pickle files with pd.DF efficiencies for a ATLAS DT recast')
    parser.add_argument('input_Dir', metavar='input_dir_path', help='Path to the Analysis directory containing cutflow files in txt') #required=True if not positional
    parser.add_argument('-o','--output', metavar='Output_dir_path', help='Path to the output directory to write the pickle files', default=".")
    args = parser.parse_args()

    inputdir = args.input_Dir
    output = args.output

    filelist = []
    results = [] # Initialize list for UL results and dictionary for Cutflow output
    masstaupairs = [(100.0, 0.012924465962305568), (100.0, 0.020483898119853464), (100.0, 0.03246479068522276), (100.0, 0.05145322575168453), (100.0, 0.08154786722400967), (100.0, 0.1292446596230557), (100.0, 0.20483898119853475), (100.0, 0.32464790685222766), (100.0, 0.5145322575168454), (100.0, 0.8154786722400966), (100.0, 1.292446596230557), (100.0, 2.048389811985348), (100.0, 3.2464790685222775), (100.0, 5.145322575168454), (100.0, 8.154786722400967), (200.0, 0.012924465962305568), (200.0, 0.020483898119853464), (200.0, 0.03246479068522276), (200.0, 0.05145322575168453), (200.0, 0.08154786722400967), (200.0, 0.1292446596230557), (200.0, 0.20483898119853475), (200.0, 0.32464790685222766), (200.0, 0.5145322575168454), (200.0, 0.8154786722400966), (200.0, 1.292446596230557), (200.0, 2.048389811985348), (200.0, 3.2464790685222775), (200.0, 5.145322575168454), (200.0, 8.154786722400967), (300.0, 0.012924465962305568), (300.0, 0.020483898119853464), (300.0, 0.03246479068522276), (300.0, 0.05145322575168453), (300.0, 0.08154786722400967), (300.0, 0.1292446596230557), (300.0, 0.20483898119853475), (300.0, 0.32464790685222766), (300.0, 0.5145322575168454), (300.0, 0.8154786722400966), (300.0, 1.292446596230557), (300.0, 2.048389811985348), (300.0, 3.2464790685222775), (300.0, 5.145322575168454), (300.0, 8.154786722400967), (400.0, 0.012924465962305568), (400.0, 0.020483898119853464), (400.0, 0.03246479068522276), (400.0, 0.05145322575168453), (400.0, 0.08154786722400967), (400.0, 0.1292446596230557), (400.0, 0.20483898119853475), (400.0, 0.32464790685222766), (400.0, 0.5145322575168454), (400.0, 0.8154786722400966), (400.0, 1.292446596230557), (400.0, 2.048389811985348), (400.0, 3.2464790685222775), (400.0, 5.145322575168454), (400.0, 8.154786722400967), (500.0, 0.012924465962305568), (500.0, 0.020483898119853464), (500.0, 0.03246479068522276), (500.0, 0.05145322575168453), (500.0, 0.08154786722400967), (500.0, 0.1292446596230557), (500.0, 0.20483898119853475), (500.0, 0.32464790685222766), (500.0, 0.5145322575168454), (500.0, 0.8154786722400966), (500.0, 1.292446596230557), (500.0, 2.048389811985348), (500.0, 3.2464790685222775), (500.0, 5.145322575168454), (500.0, 8.154786722400967), (600.0, 0.012924465962305568), (600.0, 0.020483898119853464), (600.0, 0.03246479068522276), (600.0, 0.05145322575168453), (600.0, 0.08154786722400967), (600.0, 0.1292446596230557), (600.0, 0.20483898119853475), (600.0, 0.32464790685222766), (600.0, 0.5145322575168454), (600.0, 0.8154786722400966), (600.0, 1.292446596230557), (600.0, 2.048389811985348), (600.0, 3.2464790685222775), (600.0, 5.145322575168454), (600.0, 8.154786722400967), (700.0, 0.012924465962305568), (700.0, 0.020483898119853464), (700.0, 0.03246479068522276), (700.0, 0.05145322575168453), (700.0, 0.08154786722400967), (700.0, 0.1292446596230557), (700.0, 0.20483898119853475), (700.0, 0.32464790685222766), (700.0, 0.5145322575168454), (700.0, 0.8154786722400966), (700.0, 1.292446596230557), (700.0, 2.048389811985348), (700.0, 3.2464790685222775), (700.0, 5.145322575168454), (700.0, 8.154786722400967), (800.0, 0.012924465962305568), (800.0, 0.020483898119853464), (800.0, 0.03246479068522276), (800.0, 0.05145322575168453), (800.0, 0.08154786722400967), (800.0, 0.1292446596230557), (800.0, 0.20483898119853475), (800.0, 0.32464790685222766), (800.0, 0.5145322575168454), (800.0, 0.8154786722400966), (800.0, 1.292446596230557), (800.0, 2.048389811985348), (800.0, 3.2464790685222775), (800.0, 5.145322575168454), (800.0, 8.154786722400967), (900.0, 0.012924465962305568), (900.0, 0.020483898119853464), (900.0, 0.03246479068522276), (900.0, 0.05145322575168453), (900.0, 0.08154786722400967), (900.0, 0.1292446596230557), (900.0, 0.20483898119853475), (900.0, 0.32464790685222766), (900.0, 0.5145322575168454), (900.0, 0.8154786722400966), (900.0, 1.292446596230557), (900.0, 2.048389811985348), (900.0, 3.2464790685222775), (900.0, 5.145322575168454), (900.0, 8.154786722400967), (1000.0, 0.012924465962305568), (1000.0, 0.020483898119853464), (1000.0, 0.03246479068522276), (1000.0, 0.05145322575168453), (1000.0, 0.08154786722400967), (1000.0, 0.1292446596230557), (1000.0, 0.20483898119853475), (1000.0, 0.32464790685222766), (1000.0, 0.5145322575168454), (1000.0, 0.8154786722400966), (1000.0, 1.292446596230557), (1000.0, 2.048389811985348), (1000.0, 3.2464790685222775), (1000.0, 5.145322575168454), (1000.0, 8.154786722400967)]

    #mass_pat = re.compile(r'.*wino_(\d+)GeV_.*.root')
    #tau_pat = re.compile(r'.*_(\d+\.\d)ns.*.root') #Compile regex patterns once to parse parameters from MG5 output for each point
    total_pat = re.compile(r'All : (\d+)')
    event_pat = re.compile(r'Kinematic : (\d+\.\d+)')
    final_pat = re.compile(r'0.1 < \|eta\| < 1.9 : (\d+\.\d+)')
    run_pat = re.compile(r'run_(\d+)')

    print("Arguments read: %s and %s" %(inputdir, output))
    write_path = Path(args.output) #Initialize path to folder to write results, and create any parents if needed
    write_path.mkdir(parents=True, exist_ok=True)

    #for base_path in glob.glob(inputdir+"/**/", recursive=False):
    for base_path in Path(inputdir).glob('**/*'):
        filelist.append(Path(base_path)) #
    filelist.sort(key=lambda x: runmatcher(run_pat,str(x))) #
    #print(filelist)
    for run,base_path in enumerate(filelist): #
        in_path=Path(base_path) # Initialize paths to Input and Output dirs in the MA5 results folder for the point

        print(in_path)

        mch_match=tot_match=evt_match=SR_match=tau_match=-1

        mch=masstaupairs[run][0]
        tau=masstaupairs[run][1]
        tot=[]
        evt=[]
        SR=[]

        with in_path.open('r') as infile:
            for line in infile:          #Loop over file to match the regex patters line by line, stops once all are matched for performance
                #if mch==-1:
                #    mch_match = re.search(mass_pat, line)
                if len(tot)<2:
                    tot_match = re.search(total_pat, line)
                if len(evt)<2:
                    evt_match = re.search(event_pat, line)
                if len(SR)<2:
                    SR_match = re.search(final_pat, line)
                #if tau==-1:
                #    tau_match = re.search(tau_pat, line)

                #if mch_match is not None:
                #    mch = float(mch_match[1])
                if tot_match is not None:
                    tot.append(int(tot_match[1]))
                if evt_match is not None:
                    evt.append(float(evt_match[1]))
                if SR_match is not None:
                    SR.append(float(SR_match[1]))
                #if tau_match is not None:
                #    tau = float(tau_match[1])

                #if msb!=-1 and mneu1!=-1 and ctau!=-1:
                #   break
        nor = normalizer()
        ratios = nor.ratios(mch)
        norm_events = nor.normalize_evts(mch, {'pp2chgchg':SR[0], 'pp2chgn1':SR[1]}, {'pp2chgchg':tot[0], 'pp2chgn1':tot[1]})
        results.append({'mass_ch': mch, 'tau_ns':tau, 'Total (weighted)':(w_tot:=nor.nc(mch, tot[0], tot[1])), 'Kin':(w_kin:=nor.nc(mch, evt[0], evt[1])), 'SR':(w_SR:=nor.nc(mch, SR[0],SR[1])),
                             'Evt_eff':w_kin/w_tot, 'SR_eff':w_SR/w_tot, 'Norm_evts':norm_events['pp2chgchg']+norm_events['pp2chgn1'],
                             'EvtAcc_cc':evt[0]/tot[0],'EvtAcc_cn':evt[1]/tot[1],'SReff_cc': SR[0]/tot[0], 'SReff_cn':SR[1]/tot[1]})
        #results['EW', mch, tau] = {'Total':tot[0]+tot[1], 'Kin':evt[0]+evt[1], 'SR':SR[0]+SR[1], 'Evt_eff':(evt[0]+evt[1])/(tot[0]+tot[1]), 'SR_eff':SR[0]/tot[0]} #Missing the relative cross section, need to interpolate from ref.
    pickle_path = write_path / "atlas_dt_effmap.pcl"
    with open(pickle_path,'wb') as handle:
        pickle.dump(pd.DataFrame(results), handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("Parsing has finished!")