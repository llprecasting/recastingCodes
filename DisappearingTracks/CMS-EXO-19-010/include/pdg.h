// Following contains some code from PIDCodes.h and PIDUtils.h from MCUtils (hence copyright notice below)
// with subtractions, modifications and additions by Mark Goodsell <goodsell@lpthe.jussieu.fr>
// In accordance with the licence, the namespace has been changed to HEP.

// Notable changes include the full list of particle codes that exist within pythia so that
// the charge can be looked up.
// It is not necessarily especially efficient.

// If you have a new particle not in the list below then you will need to add it for
// the HEPMC mode to work! (Especially if it is charged ... )


// -*- C++ -*-
//
// This file is part of MCUtils -- https://bitbucket.org/andybuckley/mcutils
// Copyright (C) 2013-2017 Andy Buckley <andy.buckley@cern.ch>
//
// Embedding of MCUtils code in other projects is permitted provided this
// notice is retained and the MCUtils namespace and include path are changed.
//
#pragma once
#include <algorithm>
namespace PDG {
    

/// @name Leptons
    //@{
    static const int ELECTRON = 11;
    static const int POSITRON = -ELECTRON;
    static const int EMINUS = ELECTRON;
    static const int EPLUS = POSITRON;
    static const int MUON = 13;
    static const int ANTIMUON = -MUON;
    static const int TAU = 15;
    static const int ANTITAU = -TAU;
    //@}

    /// @name Neutrinos
    //@{
    static const int NU_E = 12;
    static const int NU_EBAR = -NU_E;
    static const int NU_MU = 14;
    static const int NU_MUBAR = -NU_MU;
    static const int NU_TAU = 16;
    static const int NU_TAUBAR = -NU_TAU;
    //@}

    /// @name Bosons
    //@{
    static const int PHOTON = 22;
    static const int GAMMA = PHOTON;
    static const int GLUON = 21;
    static const int WPLUSBOSON = 24;
    static const int WMINUSBOSON = -WPLUSBOSON;
    static const int WPLUS = WPLUSBOSON;
    static const int WMINUS = WMINUSBOSON;
    static const int Z0BOSON = 23;
    static const int ZBOSON = Z0BOSON;
    static const int Z0 = Z0BOSON;
    static const int HIGGSBOSON = 25;
    static const int HIGGS = HIGGSBOSON;
    //@}

    /// @name Quarks
    //@{
    static const int DQUARK = 1;
    static const int UQUARK = 2;
    static const int SQUARK = 3;
    static const int CQUARK = 4;
    static const int BQUARK = 5;
    static const int TQUARK = 6;
    //@}

    /// @name Nucleons
    //@{
    static const int PROTON = 2212;
    static const int ANTIPROTON = -PROTON;
    static const int PBAR = ANTIPROTON;
    static const int NEUTRON = 2112;
    static const int ANTINEUTRON = -NEUTRON;
    //@}

    /// @name Light mesons
    //@{
    static const int PI0 = 111;
    static const int PIPLUS = 211;
    static const int PIMINUS = -PIPLUS;
    static const int K0L = 130;
    static const int K0S = 310;
    static const int KPLUS = 321;
    static const int KMINUS = -KPLUS;
    static const int ETA = 221;
    static const int ETAPRIME = 331;
    static const int PHI = 333;
    static const int OMEGA = 223;
    //@}

    /// @name Charmonia
    //@{
    static const int ETAC = 441;
    static const int JPSI = 443;
    static const int PSI2S = 100443;
    //@}

    /// @name Charm mesons
    //@{
    static const int D0 = 421;
    static const int DPLUS = 411;
    static const int DMINUS = -DPLUS;
    static const int DSPLUS = 431;
    static const int DSMINUS = -DSPLUS;
    //@}

    /// @name Bottomonia
    //@{
    static const int ETAB = 551;
    static const int UPSILON1S = 553;
    static const int UPSILON2S = 100553;
    static const int UPSILON3S = 200553;
    static const int UPSILON4S = 300553;
    //@}

    /// @name b mesons
    //@{
    static const int B0 = 511;
    static const int BPLUS = 521;
    static const int BMINUS = -BPLUS;
    static const int B0S = 531;
    static const int BCPLUS = 541;
    static const int BCMINUS = -BCPLUS;
    //@}

    /// @name Baryons
    //@{
    static const int LAMBDA = 3122;
    static const int SIGMA0 = 3212;
    static const int SIGMAPLUS = 3222;
    static const int SIGMAMINUS = 3112;
    static const int LAMBDACPLUS = 4122;
    static const int LAMBDACMINUS = 4122;
    static const int LAMBDAB = 5122;
    static const int XI0 = 3322;
    static const int XIMINUS = 3312;
    static const int XIPLUS = -XIMINUS;
    static const int OMEGAMINUS = 3334;
    static const int OMEGAPLUS = -OMEGAMINUS;
    //@}

    /// @name Exotic/weird stuff
    //@{
    static const int REGGEON = 110;
    static const int POMERON = 990;
    static const int ODDERON = 9990;
    static const int GRAVITON = 39;
    static const int NEUTRALINO1 = 1000022;
    static const int GRAVITINO = 1000039;
    static const int GLUINO = 1000021;
    static const int BPRIME = 7;
    static const int TPRIME = 8;
    static const int LPRIME = 17;
    static const int NUPRIME = 18;
    // static const int DARKMATTERSCALAR = 1000051;
    // static const int DARKMATTERFERMION = 1000052;
    // static const int DARKMATTERVECTOR = 1000053;
    /// @todo Add axion, black hole remnant, etc. on demand
    //@}


    /// Determine if the PID is that of a photon
    inline bool isPhoton(int pid) {
      return pid == PHOTON;
    }

    /// Determine if the PID is that of an electron or positron
    inline bool isElectron(int pid) {
      return abs(pid) == ELECTRON;
    }

    /// Determine if the PID is that of an muon or antimuon
    inline bool isMuon(int pid) {
      return abs(pid) == MUON;
    }

    /// Determine if the PID is that of an tau or antitau
    inline bool isTau(int pid) {
      return abs(pid) == TAU;
    }

    /// Determine if the PID is that of a charged lepton
    inline bool isChargedLepton(int pid) {
      const long apid = abs(pid);
      return apid == 11 || apid == 13 || apid == 15 || apid == 17;
    }
    

    /// Determine if the PID is that of a neutrino
    inline bool isNeutrino(int pid) {
      const long apid = abs(pid);
      return apid == 12 || apid == 14 || apid == 16 || apid == 18;
    }

    inline bool isLepton(int pid){
        const long apid = abs(pid);
        return (apid > 10) && (apid < 19);
    }

    inline bool isHadron(int pid){ // rather poor approximation but whatever, we assume that the user takes PDG numbers > 1000000
        const long apid = abs(pid);
        return (apid < 10) || (apid == 21) || ((apid > 100) && (apid < 1000000)) || (apid > 9000000);
    }

    inline bool isInvisible(int pid)
    {
      const long apid = abs(pid);
      return apid == 12 || apid == 14 || apid == 16 || apid == 18 || apid == 1000022 || apid == 1000023 || apid == 1000025 || apid == 1000035 ;

    }


    static const std::vector<int> neutrals={21,22,23,25,12,14,16,111,221,9000221,113,223,
    331,9010221,9000111,333,10223,10113,20113,225,20223,100221,100111,
    115,10221,9000113,9020221,10333,20333,1000223,10111,100113,100331,
    9030221,335,9010113,9020113,10225,30223,227,10115,100333,117,30113,
    9000115,10331,9010111,337,9050225,119,9060225,229,9080225,9090225,
    311,310,130,9000311,313,10313,20313,100313,10311,315,30313,10315,
    317,20315,319,421,423,10421,10423,425,511,513,515,531,533,535,441,
    443,10441,20443,10443,445,100441,100443,30443,100445,9000443,9010443,
    9020443,553,10551,20553,10553,555,100553,30553,110551,120553,100555,
    200553,300553,9000553,9010553,2112,12112,1214,22112,32112,2116,12116,
    21214,42112,31214,1218,2114,32114,1212,12114,11212,1216,21212,22114,
    11216,2118,3122,13122,3124,23122,33122,13124,43122,53122,3126,13126,
    23124,3128,23126,3212,3214,13212,13214,23212,3216,13216,23214,3218,
    3322,3324,203322,13324,103326,203326,4112,4114,4132,4312,4314,104314,
    104312,4332,4334,5122,5232};

    static const std::vector<int> chargedps={1,2,3,4,5,6,7,8,11,13,15,17,24,
      34,37,42,56,57,59,94,211,213,215,321,323,325,411,413,415,431,
      433,435,521,523,525,541,543,545,1103,1114,2101,2103,2203,2212,
      2214,2224,3101,3103,3112,3114,3201,3203,3222,3224,3303,3312,
      3314,3334,4101,4103,4122,4124,4201,4203,4212,4214,4222,4224,
      4232,4301,4303,4322,4324,4403,4412,4414,4422,4424,4432,4434,
      4444,5101,5103,5112,5114,5132,5201,5203,5222,5224,5242,5301,
      5303,5312,5314,5332,5334,5401,5403,5422,5424,5442,5444,5503,
      5512,5514,5532,5534,5554,10211,10213,10321,10323,10411,10413,
      10431,10433,10521,10523,10541,10543,14122,20213,20323,20413,
      20433,20523,20543,30323,100213,1000001,1000002,1000003,1000004,
      1000005,1000006,1000011,1000013,1000015,1000024,1000037,1000522,
      1000542,1000612,1000632,1000652,1005113,1005223,1005311,1005313,
      1005333,1006211,1006213,1006223,1006321,1006323,1009213,1009323,
      1009413,1009433,1009523,1009543,1091114,1092214,1092224,1093114,
      1093224,1093314,1093334,1094214,1094224,1094324,1095114,1095224,
      1095314,1095334,2000001,2000002,2000003,2000004,2000005,2000006,
      2000011,2000013,2000015,3000211,3000213,4000001,4000002,4000003,
      4000004,4000005,4000006,4000011,4000013,4000015,4900001,4900002,
      4900003,4900004,4900005,4900006,4900011,4900013,4900015,9000211,
      9900024,9900041,9900042,9900210,9902210};


    static const std::vector<int> charges={-1,2,-1,2,-1,2,-1,2,-3,-3,-3,-3,
      3,3,3,-1,3,3,6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,-2,-3,1,
      1,4,3,3,6,-2,-2,-3,-3,1,1,3,3,-2,-3,-3,-3,1,1,3,3,4,4,3,3,6,
      6,3,1,1,3,3,4,3,3,6,6,3,3,6,-2,-2,-3,-3,-3,1,1,3,3,3,-2,-2,-3,
      -3,-3,-3,1,1,3,3,3,3,-2,-3,-3,-3,-3,-3,3,3,3,3,3,3,3,3,3,3,3,
      3,3,3,3,3,3,3,3,3,3,-1,2,-1,2,-1,2,-3,-3,-3,3,3,-3,-3,3,3,3,
      -3,3,-3,-3,-3,3,3,6,3,3,3,3,3,3,3,3,-3,3,6,-3,3,-3,-3,3,6,3,
      -3,3,-3,-3,-1,2,-1,2,-1,2,-3,-3,-3,3,3,-1,2,-1,2,-1,2,-3,-3,
      -3,-1,2,-1,2,-1,2,-3,-3,-3,3,3,6,6,3,3}; 
 //   static std::vector<int> chargedps={24,11,13,15,1,2,3,4,5,6,211,213,9000211,10213,
 //   20213,100211,215,9000213,10211,100213,9010213,9020213,10215,217,30213,
  //  9000215,9010211,219,321,9000321,323,10323,20323,100323,10321,325,30323,
 //   10325,327,20325,329,411,413,10411,415,431,433,10431,20433,10433,435,521,
 //   523,525,541,2212,12212,2124,22212,32212,2216,12216,22124,42212,32124,2128,
 //   1114,2214,2224,31114,32214,32224,1112,2122,2222,11114,12214,12224,11112,
  //  12122,12222,1116,2126,2226,21112,22122,22222,21114,22214,22224,11116,12126,
  //  12226,1118,2218,2228,3222,3112,3114,3224,13112,13222,13114,13224,23112,
 //   23222,3116,3226,13116,13226,23114,23224,3118,3228,3312,3314,203312,13314,
 //   103316,203316,3334,203338,4122,14122,104122,204126,4212,4222,4214,4224,4232,
 //   4322,4324,104324,104322,5112,5222,5114,5224,5132,5332};

  
    inline int get3Q(int pid) {
      int apid = abs(pid);
      auto result1 = std::find(chargedps.begin(),chargedps.end(),apid);
      if(result1 == chargedps.end()) return 0;

      //int posch= charges[std::distance(chargedps,result1)];
      int index = result1-chargedps.begin();
      int posch = charges[index];
      if(pid > 0) return posch;
      return -posch;
      

    }

    inline bool isCharged(int pid){
        int apid = abs(pid);

        if(std::find(chargedps.begin(),chargedps.end(),apid) != chargedps.end())
        {
            return true;
        }
        return false;
    }

}
