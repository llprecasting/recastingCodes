// PARAMETRIZED EFFICIENCIES (https://arxiv.org/abs/1710.04901)
// https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/
// Event selection efficiency based on truth MET

#include "parametrized_truthEff.h"

namespace Pythia8{

double eventEff_MET(double Rdecay, double MET){
  double eff=1.0;
  //Before Calo
  if((0<Rdecay) && (Rdecay<1150)){
    if((200<MET) && (MET<250)){
      eff = 0.2211;
    }
    if((250<MET) && (MET<500)){
      eff = 0.8361;
    }
    if((500<MET) && (MET<2500)){
      eff = 0.9823;
    }
  }
  //Inside Calo
  if((1150<Rdecay) && (Rdecay<3870)){
    if((200<MET) && (MET<300)){
      eff = 0.4701;
    }
    if((300<MET) && (MET<1000)){
      eff = 0.7423;
    }
    if((1000<MET) && (MET<2500)){
      eff = 0.9173;
    }
  }
  //After Calo
  if(3870<Rdecay){
    if((200<MET) && (MET<250)){
      eff = 0.7351;
    }
    if((250<MET) && (MET<350)){
      eff = 0.8349;
    }
    if((350<MET) && (MET<2500)){
      eff = 0.8762;
    }
  }

  return eff;
}

// Vertex selection efficiency based on decay regions in detector
double vertexEff_Regions(double Rdecay, double massDV, int nTrk){
  double eff=1.0;
  //Region 0
  if((4.<Rdecay) && (Rdecay<22)){
    if((nTrk==5) || (nTrk==6)){
     if((10.<massDV) && (massDV<15.))  eff = 0.2802;
     if((15.<massDV) && (massDV<20.))  eff = 0.2851;
     if((20.<massDV) && (massDV<30.))  eff = 0.2774;
     if((30.<massDV) && (massDV<50.))  eff = 0.2655;
     if((50.<massDV) && (massDV<100.)) eff = 0.1954;
     if((100.<massDV)&& (massDV<200.)) eff = 0.2125;
     if((200.<massDV)&& (massDV<500.)) eff = 0.1364;
   }  
   if((nTrk==6) || (nTrk==7)){
    if((10.<massDV)  && (massDV<15.))   eff = 0.3737;
    if((15.<massDV)  && (massDV<20.))   eff = 0.4521;
    if((20.<massDV)  && (massDV<30.))   eff = 0.4458;
    if((30.<massDV)  && (massDV<50.))   eff = 0.426;
    if((50.<massDV)  && (massDV<100.))  eff = 0.3844;
    if((100.<massDV) && (massDV<200.))  eff = 0.2568;
    if((200.<massDV) && (massDV<500.))  eff = 0.2308;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2857;
   }  
   if((7<nTrk) && (nTrk<=10)){
    if((10.< massDV) && (massDV<15.))   eff = 0.4115;
    if((15.< massDV) && (massDV<20.))   eff = 0.5969;
    if((20.< massDV) && (massDV<30.))   eff = 0.6455;
    if((30.< massDV) && (massDV<50.))   eff = 0.6645;
    if((50.< massDV) && (massDV<100.))  eff = 0.6119;
    if((100.<massDV) && (massDV<200.))  eff = 0.4888;
    if((200.<massDV) && (massDV<500.))  eff = 0.4231;
    if((500.<massDV) && (massDV<5000.)) eff = 0.371;
   }  
   if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff = 0.4033;
    if((15. <massDV) && (massDV<20.))  eff = 0.6468;
    if((20. <massDV) && (massDV<30.))  eff = 0.7599;
    if((30. <massDV) && (massDV<50.))  eff = 0.8011;
    if((50. <massDV) && (massDV<100.)) eff = 0.7655;
    if((100.<massDV) && (massDV<200.)) eff = 0.7114;
    if((200.<massDV) && (massDV<500.)) eff = 0.6726;
    if((500.<massDV) && (massDV<5000.))eff = 0.6282;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff = 0.3543;
    if((15. <massDV) && (massDV<20.))  eff = 0.6216;
    if((20. <massDV) && (massDV<30.))  eff = 0.7561;
    if((30. <massDV) && (massDV<50.))  eff = 0.8215;
    if((50. <massDV) && (massDV<100.)) eff = 0.8217;
    if((100.<massDV) && (massDV<200.)) eff = 0.8138;
    if((200.<massDV) && (massDV<500.)) eff = 0.791;
    if((500.<massDV) && (massDV<5000.))eff = 0.7518;
   } 
   if((20<nTrk) && (nTrk<=30)){
    if((10.< massDV) && (massDV<15.))  eff = 0.2766;
    if((15.< massDV) && (massDV<20.))  eff = 0.5018;
    if((20.< massDV) && (massDV<30.))  eff = 0.678;
    if((30.< massDV) && (massDV<50.))  eff = 0.809;
    if((50.< massDV) && (massDV<100.)) eff = 0.845;
    if((100.<massDV) && (massDV<200.)) eff = 0.8591;
    if((200.<massDV) && (massDV<500.)) eff = 0.8581;
    if((500.<massDV) && (massDV<5000.)) eff = 0.8383;
   }
   if((30<nTrk) && (nTrk<=50)){
    if((20.< massDV) && (massDV<30.))  eff = 0.5882;
    if((30.< massDV) && (massDV<50.))  eff = 0.7438;
    if((50.< massDV) && (massDV<100.)) eff = 0.87;
    if((100.<massDV) && (massDV<200.)) eff = 0.8823;
    if((200.<massDV) && (massDV<500.)) eff = 0.8855;
    if((500.<massDV) && (massDV<5000.)) eff = 0.8836;
   }     
   if((50<nTrk) && (nTrk<=200)){
    if((30. <massDV) && (massDV<50.))  eff = 1.;
    if((50. <massDV) && (massDV<100.)) eff = 0.84;
    if((100.<massDV) && (massDV<200.)) eff = 0.89;
    if((200.<massDV) && (massDV<500.)) eff = 0.8928;
    if((500.<massDV) && (massDV<5000.)) eff = 0.8934;
   }  
  }
  //Region 1
  if((22.<Rdecay) && (Rdecay<25.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.< massDV) && (massDV<15.)) eff = 0.06298;
      if((15.< massDV) && (massDV<20.)) eff = 0.07806;
      if((20.< massDV) && (massDV<30.)) eff = 0.07449;
      if((30.< massDV) && (massDV<50.)) eff = 0.04;
      if((50.< massDV) && (massDV<100.)) eff = 0.01613;
      if((100.<massDV) && (massDV<200.)) eff = 0.125;
   }  
   if((nTrk==6) || (nTrk==7)){
      if((10. <massDV) && (massDV<15.)) eff = 0.1039;
      if((15. <massDV) && (massDV<20.)) eff = 0.1169;
      if((20. <massDV) && (massDV<30.)) eff = 0.1073;
      if((30. <massDV) && (massDV<50.)) eff = 0.09719;
      if((50. <massDV) && (massDV<100.)) eff = 0.07746;
      if((100.<massDV) && (massDV<200.)) eff = 0.1111;
   }  
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.1128;
    if((15. <massDV) && (massDV<20.)) eff = 0.1577;
    if((20. <massDV) && (massDV<30.)) eff = 0.1887;
    if((30. <massDV) && (massDV<50.)) eff = 0.1834;
    if((50. <massDV) && (massDV<100.)) eff = 0.1533;
    if((100.<massDV) && (massDV<200.)) eff = 0.1278;
    if((200.<massDV) && (massDV<500.)) eff = 0.129;
   }  
   if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff = 0.1085;
    if((15. <massDV) && (massDV<20.))  eff = 0.1812;
    if((20. <massDV) && (massDV<30.))  eff = 0.2364;
    if((30. <massDV) && (massDV<50.))  eff = 0.246;
    if((50. <massDV) && (massDV<100.)) eff = 0.2322;
    if((100.<massDV) && (massDV<200.)) eff = 0.2019;
    if((200.<massDV) && (massDV<500.)) eff = 0.195;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1299;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff = 0.1463;
    if((15. <massDV) && (massDV<20.))  eff = 0.1989;
    if((20. <massDV) && (massDV<30.))  eff = 0.2282;
    if((30. <massDV) && (massDV<50.))  eff = 0.2386;
    if((50. <massDV) && (massDV<100.)) eff = 0.2558;
    if((100.<massDV) && (massDV<200.)) eff = 0.2557;
    if((200.<massDV) && (massDV<500.)) eff = 0.2399;
    if((500.<massDV) && (massDV<5000.)) eff = 0.205;
   }   
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.)) eff = 0.1053;
    if((20. <massDV) && (massDV<30.)) eff = 0.1442;
    if((30. <massDV) && (massDV<50.)) eff = 0.2381;
    if((50. <massDV) && (massDV<100.)) eff = 0.2723;
    if((100.<massDV) && (massDV<200.)) eff = 0.2785;
    if((200.<massDV) && (massDV<500.)) eff = 0.2776;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2335;
   }   
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff = 0.1111;
    if((50. <massDV) && (massDV<100.)) eff = 0.2438;
    if((100.<massDV) && (massDV<200.)) eff = 0.2757;
    if((200.<massDV) && (massDV<500.)) eff = 0.2942;
    if((500.<massDV) && (massDV<5000.)) eff = 0.3022;
  }
   if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV) && (massDV<100.)) eff = 0.5;
    if((100.<massDV) && (massDV<200.)) eff = 0.2692;
    if((200.<massDV) && (massDV<500.)) eff = 0.3011;
    if((500.<massDV) && (massDV<5000.)) eff = 0.29;
   }
  }
  //Region 2
  if((25.<Rdecay) && (Rdecay<29.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.< massDV) && (massDV<15.)) eff = 0.08709;
      if((15.< massDV) && (massDV<20.)) eff = 0.09924;
      if((20.< massDV) && (massDV<30.)) eff = 0.07775;
      if((30.< massDV) && (massDV<50.)) eff = 0.1101;
      if((50.< massDV) && (massDV<100.)) eff = 0.1585;
      if((100.< massDV) && (massDV<200.)) eff = 0.1111;
    }  
    if((nTrk==6) || (nTrk==7)){
      if((10. < massDV) && (massDV<15.)) eff = 0.1407;
      if((15. < massDV) && (massDV<20.)) eff = 0.1804;
      if((20. < massDV) && (massDV<30.)) eff = 0.2052;
      if((30. < massDV) && (massDV<50.)) eff = 0.1549;
      if((50. < massDV) && (massDV<100.)) eff = 0.1504;
      if((100.< massDV) && (massDV<200.)) eff = 0.125;
      if((200.< massDV) && (massDV<500.)) eff = 0.25;
    } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff   =0.1602;
    if((15. <massDV) && (massDV<20.)) eff   =0.2565;
    if((20. <massDV) && (massDV<30.)) eff   =0.2773;
    if((30. <massDV) && (massDV<50.)) eff   =0.3059;
    if((50. <massDV) && (massDV<100.)) eff  =0.2537;
    if((100.<massDV) && (massDV<200.)) eff =0.2091;
    if((200.<massDV) && (massDV<500.)) eff =0.1711;
    }  
   if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff = 0.1699;
    if((15. <massDV) && (massDV<20.))  eff = 0.2704;
    if((20. <massDV) && (massDV<30.))  eff = 0.3847;
    if((30. <massDV) && (massDV<50.))  eff = 0.4098;
    if((50. <massDV) && (massDV<100.)) eff = 0.3943;
    if((100.<massDV) && (massDV<200.)) eff = 0.3351;
    if((200.<massDV) && (massDV<500.)) eff = 0.3103;
    if((500.<massDV) && (massDV<5000.)) eff =0.1932;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff = 0.08537;
    if((15. <massDV) && (massDV<20.))  eff = 0.295;
    if((20. <massDV) && (massDV<30.))  eff = 0.4156;
    if((30. <massDV) && (massDV<50.))  eff = 0.4211;
    if((50. <massDV) && (massDV<100.)) eff =  0.4452;
    if((100.<massDV) && (massDV<200.)) eff = 0.4396;
    if((200.<massDV) && (massDV<500.)) eff =  0.3898;
    if((500.<massDV) && (massDV<5000.)) eff = 0.3735;
    }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.))  eff =  0.125;
    if((20. <massDV) && (massDV<30.))  eff =  0.3211;
    if((30. <massDV) && (massDV<50.))  eff =  0.3761;
    if((50. <massDV) && (massDV<100.)) eff =  0.4242;
    if((100.<massDV) && (massDV<200.)) eff =   0.4816;
    if((200.<massDV) && (massDV<500.)) eff =   0.477;
    if((500.<massDV) && (massDV<5000.)) eff =  0.4318;
   }
  if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff =  0.3077;
    if((50. <massDV) && (massDV<100.)) eff =  0.449;
    if((100.<massDV) && (massDV<200.)) eff =  0.4882;
    if((200.<massDV) && (massDV<500.)) eff =  0.5007;
    if((500.<massDV) && (massDV<5000.)) eff = 0.4932;
   }
  if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV) && (massDV<100.)) eff =  1.0;
    if((100.<massDV) && (massDV<200.)) eff =  0.4673;
    if((200.<massDV) && (massDV<500.)) eff =  0.5015;
    if((500.<massDV) && (massDV<5000.)) eff = 0.4849;
   }
  }
 //Region 3
  if((29.<Rdecay) && (Rdecay<38.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =   0.01048;
      if((15.<massDV) && (massDV<20.)) eff =   0.006224;
      if((20.<massDV) && (massDV<30.)) eff =   0.01079;
      if((30.<massDV) && (massDV<50.)) eff =   0.01144;
      if((50.<massDV) && (massDV<100.)) eff =  0.01399;
   } 
    if((nTrk==6) || (nTrk==7)){
      if((10.<massDV) && (massDV<15.)) eff =  0.01149;
      if((15.<massDV) && (massDV<20.)) eff =  0.01918;
      if((20.<massDV) && (massDV<30.)) eff =  0.02122;
      if((30.<massDV) && (massDV<50.)) eff =  0.01764;
      if((50.<massDV) && (massDV<100.)) eff = 0.01205;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.02023;
    if((15. <massDV) && (massDV<20.)) eff = 0.02106;
    if((20. <massDV) && (massDV<30.)) eff = 0.03253;
    if((30. <massDV) && (massDV<50.)) eff = 0.03222;
    if((50. <massDV) && (massDV<100.)) eff = 0.02877;
    if((100.<massDV) && (massDV<200.)) eff = 0.01479;
    if((200.<massDV) && (massDV<500.)) eff = 0.02752;
 }
   if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff =  0.01583;
    if((15. <massDV) && (massDV<20.))  eff =  0.03566;
    if((20. <massDV) && (massDV<30.))  eff =  0.03645;
    if((30. <massDV) && (massDV<50.))  eff =  0.04561;
    if((50. <massDV) && (massDV<100.)) eff =  0.0487;
    if((100.<massDV) && (massDV<200.)) eff =  0.0422;
    if((200.<massDV) && (massDV<500.)) eff =  0.02811;
    if((500.<massDV) && (massDV<5000.)) eff = 0.04762;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff =  0.0125;
    if((15. <massDV) && (massDV<20.))  eff =  0.0142;
    if((20. <massDV) && (massDV<30.))  eff =  0.03707;
    if((30. <massDV) && (massDV<50.))  eff =  0.0464;
    if((50. <massDV) && (massDV<100.)) eff =  0.05051;
    if((100.<massDV) && (massDV<200.)) eff =  0.04665;
    if((200.<massDV) && (massDV<500.)) eff =  0.04462;
    if((500.<massDV) && (massDV<5000.)) eff = 0.03311;
 }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.))  eff =  0.02857;
    if((20. <massDV) && (massDV<30.))  eff =  0.02679;
    if((30. <massDV) && (massDV<50.))  eff =  0.03477;
    if((50. <massDV) && (massDV<100.)) eff =  0.04878;
    if((100.<massDV) && (massDV<200.)) eff =  0.05055;
    if((200.<massDV) && (massDV<500.)) eff =  0.05072;
    if((500.<massDV) && (massDV<5000.)) eff = 0.04425;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff =  0.01449;
    if((50. <massDV) && (massDV<100.)) eff =  0.05893;
    if((100.<massDV) && (massDV<200.)) eff =  0.04793;
    if((200.<massDV) && (massDV<500.)) eff =  0.05554;
    if((500.<massDV) && (massDV<5000.)) eff = 0.05854;
   }
    if((50<nTrk) && (nTrk<=200)){
    if((100.<massDV) && (massDV<200.)) eff =  0.05164;
    if((200.<massDV) && (massDV<500.)) eff =  0.05875;
    if((500.<massDV) && (massDV<5000.)) eff = 0.05403;
  }

}
 //Region 4
  if((38.<Rdecay) && (Rdecay<46.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =  0.04526;
      if((15.<massDV) && (massDV<20.)) eff =  0.05573;
      if((20.<massDV) && (massDV<30.)) eff =  0.04129;
      if((30.<massDV) && (massDV<50.)) eff =  0.02875;
      if((50.<massDV) && (massDV<100.)) eff = 0.05882;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10. <massDV) && (massDV<15.)) eff =  0.07339;
      if((15. <massDV) && (massDV<20.)) eff =  0.06595;
      if((20. <massDV) && (massDV<30.)) eff =  0.07368;
      if((30. <massDV) && (massDV<50.)) eff =  0.07339;
      if((50. <massDV) && (massDV<100.)) eff = 0.05952;
      if((100.<massDV) && (massDV<200.)) eff = 0.07692;
      if((200.<massDV) && (massDV<500.)) eff = 0.125;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff =  0.08318;
    if((15. <massDV) && (massDV<20.)) eff =  0.1195;
    if((20. <massDV) && (massDV<30.)) eff =  0.1278;
    if((30. <massDV) && (massDV<50.)) eff =  0.1332;
    if((50. <massDV) && (massDV<100.)) eff = 0.1103;
    if((100.<massDV) && (massDV<200.)) eff = 0.068;
    if((200.<massDV) && (massDV<500.)) eff = 0.04598;
 }

  if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff =  0.08327;
    if((15. <massDV) && (massDV<20.))  eff =  0.1588;
    if((20. <massDV) && (massDV<30.))  eff =  0.2098;
    if((30. <massDV) && (massDV<50.))  eff =  0.2214;
    if((50. <massDV) && (massDV<100.)) eff =  0.1871;
    if((100.<massDV) && (massDV<200.)) eff =  0.1449;
    if((200.<massDV) && (massDV<500.)) eff =  0.1491;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1339;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff = 0.09449;
    if((15. <massDV) && (massDV<20.))  eff = 0.146;
    if((20. <massDV) && (massDV<30.))  eff = 0.2238;
    if((30. <massDV) && (massDV<50.))  eff = 0.2458;
    if((50. <massDV) && (massDV<100.)) eff = 0.2425;
    if((100.<massDV) && (massDV<200.)) eff = 0.2282;
    if((200.<massDV) && (massDV<500.)) eff = 0.2031;
    if((500.<massDV) && (massDV<5000.)) eff =0.1425;
 }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.))  eff =  0.04167;
    if((20. <massDV) && (massDV<30.))  eff =  0.1302;
    if((30. <massDV) && (massDV<50.))  eff =  0.199;
    if((50. <massDV) && (massDV<100.)) eff =  0.2701;
    if((100.<massDV) && (massDV<200.)) eff =  0.2834;
    if((200.<massDV) && (massDV<500.)) eff =  0.2658;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2325;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV<50.))  eff =  0.1481;
    if((50. <massDV<100.)) eff =  0.2636;
    if((100.<massDV<200.)) eff =  0.3228;
    if((200.<massDV<500.)) eff =  0.3356;
    if((500.<massDV<5000.)) eff = 0.3272;
 }
    if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV<100.)) eff =  0.25;
    if((100.<massDV<200.)) eff =  0.3054;
    if((200.<massDV<500.)) eff =  0.3703;
    if((500.<massDV<5000.)) eff = 0.377;
  }
 }
 //Region 5
  if((46.<Rdecay) && (Rdecay<73.)){
    if((nTrk==5) || (nTrk==6)){
      if((10. <massDV) && (massDV<15.)) eff =  0.04293;
      if((15. <massDV) && (massDV<20.)) eff =  0.04807;
      if((20. <massDV) && (massDV<30.)) eff =  0.04606;
      if((30. <massDV) && (massDV<50.)) eff =  0.0498;
      if((50. <massDV) && (massDV<100.)) eff = 0.02996;
      if((200.<massDV) && (massDV<500.)) eff = 0.1429;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10. <massDV) && (massDV<15.)) eff =  0.05983;
      if((15. <massDV) && (massDV<20.)) eff =  0.07004;
      if((20. <massDV) && (massDV<30.)) eff =  0.07386;
      if((30. <massDV) && (massDV<50.)) eff =  0.0672;
      if((50. <massDV) && (massDV<100.)) eff = 0.06068;
      if((100.<massDV) && (massDV<200.)) eff = 0.01695;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff =  0.07426;
    if((15. <massDV) && (massDV<20.)) eff =  0.1059;
    if((20. <massDV) && (massDV<30.)) eff =  0.1198;
    if((30. <massDV) && (massDV<50.)) eff =  0.1189;
    if((50. <massDV) && (massDV<100.)) eff = 0.09594;
    if((100.<massDV) && (massDV<200.)) eff = 0.05743;
    if((200.<massDV) && (massDV<500.)) eff = 0.05677;
 }
  if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff =  0.087;
    if((15. <massDV) && (massDV<20.))  eff =  0.1356;
    if((20. <massDV) && (massDV<30.))  eff =  0.1665;
    if((30. <massDV) && (massDV<50.))  eff =  0.182;
    if((50. <massDV) && (massDV<100.)) eff =  0.161;
    if((100.<massDV) && (massDV<200.)) eff =  0.131;
    if((200.<massDV) && (massDV<500.)) eff =  0.1185;
    if((500.<massDV) && (massDV<5000.)) eff = 0.09375;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.))  eff = 0.04104;
    if((15. <massDV) && (massDV<20.))  eff = 0.1332;
    if((20. <massDV) && (massDV<30.))  eff = 0.1943;
    if((30. <massDV) && (massDV<50.))  eff = 0.221;
    if((50. <massDV) && (massDV<100.)) eff = 0.1947;
    if((100.<massDV) && (massDV<200.)) eff = 0.1795;
    if((200.<massDV) && (massDV<500.)) eff = 0.1493;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1384;
 }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.)) eff =  0.09459;
    if((20. <massDV) && (massDV<30.)) eff =  0.1207;
    if((30. <massDV) && (massDV<50.)) eff =  0.2087;
    if((50. <massDV) && (massDV<100.)) eff = 0.2218;
    if((100.<massDV) && (massDV<200.)) eff =  0.234;
    if((200.<massDV) && (massDV<500.)) eff =  0.2208;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1939;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff =  0.1382;
    if((50. <massDV) && (massDV<100.)) eff =  0.2279;
    if((100.<massDV) && (massDV<200.)) eff =  0.2789;
    if((200.<massDV) && (massDV<500.)) eff =  0.287;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2653;
 }
  if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV) && (massDV<100.)) eff =  0.07692;
    if((100.<massDV) && (massDV<200.)) eff =  0.2927;
    if((200.<massDV) && (massDV<500.)) eff =  0.3211;
    if((500.<massDV) && (massDV<5000.)) eff = 0.3204;
  }
 }

 //Region 6
  if((73.<Rdecay) && (Rdecay<84.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =  0.08585;
      if((15.<massDV) && (massDV<20.)) eff =  0.1027;
      if((20.<massDV) && (massDV<30.)) eff =  0.1188;
      if((30.<massDV) && (massDV<50.)) eff =  0.09314;
      if((50.<massDV) && (massDV<100.)) eff = 0.1772;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10. <massDV) && (massDV<15.)) eff = 0.141;
      if((15. <massDV) && (massDV<20.)) eff = 0.1692;
      if((20. <massDV) && (massDV<30.)) eff = 0.1656;
      if((30. <massDV) && (massDV<50.)) eff = 0.1588;
      if((50. <massDV) && (massDV<100.)) eff = 0.1719;
      if((100.<massDV) && (massDV<200.)) eff = 0.125;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.))  eff = 0.1774;
    if((15. <massDV) && (massDV<20.))  eff = 0.225;
    if((20. <massDV) && (massDV<30.))  eff = 0.2749;
    if((30. <massDV) && (massDV<50.))  eff = 0.2751;
    if((50. <massDV) && (massDV<100.)) eff = 0.2597;
    if((100.<massDV) && (massDV<200.)) eff = 0.2229;
    if((200.<massDV) && (massDV<500.)) eff = 0.08065;
    if((500.<massDV) && (massDV<5000.)) eff = 0.3;
 }
  if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.)) eff =  0.1919;
    if((15. <massDV) && (massDV<20.)) eff =  0.311;
    if((20. <massDV) && (massDV<30.)) eff =  0.3974;
    if((30. <massDV) && (massDV<50.)) eff =  0.433;
    if((50. <massDV) && (massDV<100.)) eff = 0.3883;
    if((100.<massDV) && (massDV<200.)) eff = 0.3116;
    if((200.<massDV) && (massDV<500.)) eff = 0.2369;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.)) eff =  0.1852;
    if((15. <massDV) && (massDV<20.)) eff =  0.2565;
    if((20. <massDV) && (massDV<30.)) eff =  0.4545;
    if((30. <massDV) && (massDV<50.)) eff =  0.5376;
    if((50. <massDV) && (massDV<100.)) eff = 0.4785;
    if((100.<massDV) && (massDV<200.)) eff =  0.4227;
    if((200.<massDV) && (massDV<500.)) eff =  0.3787;
    if((500.<massDV) && (massDV<5000.)) eff = 0.3006;
 }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.)) eff = 0.2353;
    if((20. <massDV) && (massDV<30.)) eff = 0.305;
    if((30. <massDV) && (massDV<50.)) eff = 0.4257;
    if((50. <massDV) && (massDV<100.)) eff = 0.5127;
    if((100.<massDV) && (massDV<200.)) eff = 0.5512;
    if((200.<massDV) && (massDV<500.)) eff = 0.5131;
    if((500.<massDV) && (massDV<5000.)) eff = 0.4621;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff =  0.2564;
    if((50. <massDV) && (massDV<100.)) eff =  0.4933;
    if((100.<massDV) && (massDV<200.)) eff =  0.6646;
    if((200.<massDV) && (massDV<500.)) eff =  0.6902;
    if((500.<massDV) && (massDV<5000.)) eff = 0.634;
 }
  if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV) && (massDV<100.)) eff =  0.5;
    if((100.<massDV) && (massDV<200.)) eff =  0.7143;
    if((200.<massDV) && (massDV<500.)) eff =  0.8095;
    if((500.<massDV) && (massDV<5000.)) eff = 0.7987;
  }
 }
 //Region 7
  if((84.<Rdecay) && (Rdecay<111.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =  0.02209;
      if((15.<massDV) && (massDV<20.)) eff =  0.0227;
      if((20.<massDV) && (massDV<30.)) eff =  0.03417;
      if((30.<massDV) && (massDV<50.)) eff =  0.02653;
      if((50.<massDV) && (massDV<100.)) eff = 0.06122;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10.<massDV) && (massDV<15.)) eff = 0.03759;
      if((15.<massDV) && (massDV<20.)) eff = 0.03088;
      if((20.<massDV) && (massDV<30.)) eff = 0.03767;
      if((30.<massDV) && (massDV<50.)) eff = 0.02207;
      if((50.<massDV) && (massDV<100.)) eff = 0.04331;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.04217;
    if((15. <massDV) && (massDV<20.)) eff = 0.06129;
    if((20. <massDV) && (massDV<30.)) eff = 0.06257;
    if((30. <massDV) && (massDV<50.)) eff = 0.06395;
    if((50. <massDV) && (massDV<100.)) eff = 0.06502;
    if((100.<massDV) && (massDV<200.)) eff = 0.01994;
    if((200.<massDV) && (massDV<500.)) eff = 0.024;
 }
  if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.))  eff = 0.03908;
    if((15. <massDV) && (massDV<20.))  eff = 0.08071;
    if((20. <massDV) && (massDV<30.))  eff = 0.1045;
    if((30. <massDV) && (massDV<50.))  eff = 0.1235;
    if((50. <massDV) && (massDV<100.)) eff = 0.1063;
    if((100.<massDV) && (massDV<200.)) eff =  0.06798;
    if((200.<massDV) && (massDV<500.)) eff =  0.05522;
    if((500.<massDV) && (massDV<5000.)) eff = 0.005263;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.)) eff = 0.0303;
    if((15. <massDV) && (massDV<20.)) eff = 0.04902;
    if((20. <massDV) && (massDV<30.)) eff = 0.1097;
    if((30. <massDV) && (massDV<50.)) eff = 0.1516;
    if((50. <massDV) && (massDV<100.)) eff = 0.1374;
    if((100.<massDV) && (massDV<200.)) eff = 0.1048;
    if((200.<massDV) && (massDV<500.)) eff = 0.09298;
    if((500.<massDV) && (massDV<5000.)) eff = 0.05498;
 }
   if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.)) eff = 0.02632;
    if((20. <massDV) && (massDV<30.)) eff = 0.075;
    if((30. <massDV) && (massDV<50.)) eff = 0.1139;
    if((50. <massDV) && (massDV<100.)) eff = 0.1593;
    if((100.<massDV) && (massDV<200.)) eff = 0.1582;
    if((200.<massDV) && (massDV<500.)) eff = 0.1416;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1159;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.))  eff = 0.01515;
    if((50. <massDV) && (massDV<100.)) eff = 0.1503;
    if((100.<massDV) && (massDV<200.)) eff = 0.2051;
    if((200.<massDV) && (massDV<500.)) eff = 0.2178;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1802;
 }
  if((50<nTrk) && (nTrk<=200)){
    if((100.<massDV) && (massDV<200.)) eff =  0.225;
    if((200.<massDV) && (massDV<500.)) eff =  0.278;
    if((500.<massDV) && (massDV<5000.)) eff = 0.270;
  }

}
//Region 8
  if((111.<Rdecay) && (Rdecay<120.)){
    if((nTrk==5) || (nTrk==6)){
      if((10. <massDV) && (massDV<15.)) eff = 0.0399;
      if((15. <massDV) && (massDV<20.)) eff = 0.02765;
      if((20. <massDV) && (massDV<30.)) eff = 0.03723;
      if((30. <massDV) && (massDV<50.)) eff = 0.05405;
      if((50. <massDV) && (massDV<100.)) eff = 0.05882;
      if((100.<massDV) && (massDV<200.)) eff = 0.1667;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10.<massDV) && (massDV<15.)) eff =  0.0532;
      if((15.<massDV) && (massDV<20.)) eff =  0.0905;
      if((20.<massDV) && (massDV<30.)) eff =  0.0606;
      if((30.<massDV) && (massDV<50.)) eff =  0.0529;
      if((50.<massDV) && (massDV<100.)) eff = 0.0684;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.0501;
    if((15. <massDV) && (massDV<20.)) eff = 0.0783;
    if((20. <massDV) && (massDV<30.)) eff = 0.09736;
    if((30. <massDV) && (massDV<50.)) eff = 0.07677;
    if((50. <massDV) && (massDV<100.)) eff = 0.08451;
    if((100.<massDV) && (massDV<200.)) eff = 0.08491;
    if((200.<massDV) && (massDV<500.)) eff = 0.07692; 
 }
   if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.)) eff = 0.04339;
    if((15. <massDV) && (massDV<20.)) eff = 0.09733;
    if((20. <massDV) && (massDV<30.)) eff = 0.1602;
    if((30. <massDV) && (massDV<50.)) eff = 0.1664;
    if((50. <massDV) && (massDV<100.)) eff = 0.1506;
    if((100.<massDV) && (massDV<200.)) eff = 0.1033;
    if((200.<massDV) && (massDV<500.)) eff = 0.05618;
    if((500.<massDV) && (massDV<5000.)) eff = 0.07843;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((15. <massDV) && (massDV<20.))  eff = 0.07273;
    if((20. <massDV) && (massDV<30.))  eff = 0.1429;
    if((30. <massDV) && (massDV<50.))  eff = 0.2091;
    if((50. <massDV) && (massDV<100.)) eff = 0.2126;
    if((100.<massDV) && (massDV<200.)) eff = 0.1778;
    if((200.<massDV) && (massDV<500.)) eff = 0.1327;
    if((500.<massDV) && (massDV<5000.)) eff = 0.07729;
 }
    if((20<nTrk) && (nTrk<=30)){
    if((20. <massDV) && (massDV<30.)) eff = 0.07895;
    if((30. <massDV) && (massDV<50.)) eff = 0.1193;
    if((50. <massDV) && (massDV<100.)) eff = 0.2051;
    if((100.<massDV) && (massDV<200.)) eff = 0.2299;
    if((200.<massDV) && (massDV<500.)) eff = 0.1977;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1584;
 }
   if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.)) eff =  0.07407;
    if((50. <massDV) && (massDV<100.)) eff = 0.1886;
    if((100.<massDV) && (massDV<200.)) eff =  0.3106;
    if((200.<massDV) && (massDV<500.)) eff =  0.3348;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2758;
 }
   if((50<nTrk) && (nTrk<=200)){
    if((50. <massDV) && (massDV<100.)) eff = 0.5;
    if((100.<massDV) && (massDV<200.)) eff = 0.2658;
    if((200.<massDV) && (massDV<500.)) eff = 0.463;
    if((500.<massDV) && (massDV<5000.)) eff = 0.4624;
  }
}
//Region 9
  if((120.<Rdecay) && (Rdecay<145.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =  0.01618;
      if((15.<massDV) && (massDV<20.)) eff =  0.01306;
      if((20.<massDV) && (massDV<30.)) eff =  0.01833;
      if((30.<massDV) && (massDV<50.)) eff =  0.01293;
      if((50.<massDV) && (massDV<100.)) eff = 0.04878;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10.<massDV) && (massDV<15.)) eff =  0.02299;
      if((15.<massDV) && (massDV<20.)) eff =  0.0379;
      if((20.<massDV) && (massDV<30.)) eff =  0.03851;
      if((30.<massDV) && (massDV<50.)) eff =  0.03326;
      if((50.<massDV) && (massDV<100.)) eff = 0.03311;
   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.02984;
    if((15. <massDV) && (massDV<20.)) eff = 0.04244;
    if((20. <massDV) && (massDV<30.)) eff = 0.05501;
    if((30. <massDV) && (massDV<50.)) eff = 0.04615;
    if((50. <massDV) && (massDV<100.)) eff =0.05382;
    if((100.<massDV) && (massDV<200.)) eff = 0.02551;
 }
    if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.)) eff = 0.02929;
    if((15. <massDV) && (massDV<20.)) eff = 0.05095;
    if((20. <massDV) && (massDV<30.)) eff = 0.07778;
    if((30. <massDV) && (massDV<50.)) eff = 0.09255;
    if((50. <massDV) && (massDV<100.)) eff =  0.07854;
    if((100.<massDV) && (massDV<200.)) eff =  0.05007;
    if((200.<massDV) && (massDV<500.)) eff =  0.03302;
    if((500.<massDV) && (massDV<5000.)) eff = 0.01869;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.)) eff = 0.01527;
    if((15. <massDV) && (massDV<20.)) eff = 0.02734;
    if((20. <massDV) && (massDV<30.)) eff = 0.08413;
    if((30. <massDV) && (massDV<50.)) eff = 0.1242;
    if((50. <massDV) && (massDV<100.)) eff = 0.1095;
    if((100.<massDV) && (massDV<200.)) eff = 0.08245;
    if((200.<massDV) && (massDV<500.)) eff = 0.06969;
    if((500.<massDV) && (massDV<5000.)) eff = 0.05581;
 }
  if((20<nTrk) && (nTrk<=30)){
    if((20. <massDV<30.)) eff = 0.04651;
    if((30. <massDV<50.)) eff = 0.1159;
    if((50. <massDV<100.)) eff = 0.1213;
    if((100.<massDV<200.)) eff = 0.1332;
    if((200.<massDV<500.)) eff = 0.107;
    if((500.<massDV<5000.)) eff = 0.06965;
 }
    if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.)) eff = 0.06383;
    if((50. <massDV) && (massDV<100.)) eff = 0.1259;
    if((100.<massDV) && (massDV<200.)) eff = 0.179;
    if((200.<massDV) && (massDV<500.)) eff = 0.1752;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1454;
 }
    if((50<nTrk) && (nTrk<=200)){
    if((100.<massDV) && (massDV<200.)) eff =  0.1827;
    if((200.<massDV) && (massDV<500.)) eff =  0.2426;
    if((500.<massDV) && (massDV<5000.)) eff = 0.2236;
  }
 }
//Region 10
  if((145.<Rdecay) && (Rdecay<180.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff = 0.05729;
      if((15.<massDV) && (massDV<20.)) eff = 0.08098;
      if((20.<massDV) && (massDV<30.)) eff = 0.06491;
      if((30.<massDV) && (massDV<50.)) eff = 0.08754;
      if((50.<massDV) && (massDV<100.)) eff = 0.1;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10. <massDV<15.)) eff = 0.0903;
      if((15. <massDV<20.)) eff = 0.1148;
      if((20. <massDV<30.)) eff = 0.09257;
      if((30. <massDV<50.)) eff = 0.1162;
      if((50. <massDV<100.)) eff = 0.15;
      if((100.<massDV<200.)) eff = 0.1364;

   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.09463;
    if((15. <massDV) && (massDV<20.)) eff = 0.1384;
    if((20. <massDV) && (massDV<30.)) eff = 0.1577;
    if((30. <massDV) && (massDV<50.)) eff = 0.171;
    if((50. <massDV) && (massDV<100.)) eff =  0.1402;
    if((100.<massDV) && (massDV<200.)) eff = 0.1056;
    if((200.<massDV) && (massDV<500.)) eff = 0.08 ;

 }
    if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.)) eff = 0.06274;
    if((15. <massDV) && (massDV<20.)) eff = 0.1812;
    if((20. <massDV) && (massDV<30.)) eff = 0.2492;
    if((30. <massDV) && (massDV<50.)) eff = 0.2839;
    if((50. <massDV) && (massDV<100.)) eff = 0.2449;
    if((100.<massDV) && (massDV<200.)) eff = 0.1481;
    if((200.<massDV) && (massDV<500.)) eff = 0.0967;
    if((500.<massDV) && (massDV<5000.)) eff = 0.06316;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.)) eff = 0.07752;
    if((15. <massDV) && (massDV<20.)) eff = 0.1029;
    if((20. <massDV) && (massDV<30.)) eff = 0.2402;
    if((30. <massDV) && (massDV<50.)) eff = 0.3456;
    if((50. <massDV) && (massDV<100.)) eff = 0.3193;
    if((100.<massDV) && (massDV<200.)) eff = 0.2494;
    if((200.<massDV) && (massDV<500.)) eff = 0.1976;
    if((500.<massDV) && (massDV<5000.)) eff = 0.1374;
 }
  if((20<nTrk) && (nTrk<=30)){
    if((15. <massDV) && (massDV<20.))  eff = 0.03448;
    if((20. <massDV) && (massDV<30.))  eff = 0.1397;
    if((30. <massDV) && (massDV<50.))  eff = 0.2534;
    if((50. <massDV) && (massDV<100.)) eff = 0.3387;
    if((100.<massDV) && (massDV<200.)) eff = 0.3675;
    if((200.<massDV) && (massDV<500.)) eff = 0.3175;
    if((500.<massDV) && (massDV<5000.)) eff = 0.254;
 }
    if((30<nTrk) && (nTrk<=50)){
    if((30. <massDV) && (massDV<50.)) eff = 0.08772;
    if((50. <massDV) && (massDV<100.)) eff = 0.2887;
    if((100.<massDV) && (massDV<200.)) eff = 0.4824;
    if((200.<massDV) && (massDV<500.)) eff = 0.5056;
    if((500.<massDV) && (massDV<5000.)) eff = 0.4423;
 }
    if((50<nTrk) && (nTrk<=200)){
    if((100.<massDV) && (massDV<200.)) eff = 0.4823;
    if((200.<massDV) && (massDV<500.)) eff = 0.6816;
    if((500.<massDV) && (massDV<5000.)) eff = 0.6378;
  }
}

//Region 11
  if((180.<Rdecay) && (Rdecay<300.)){
    if((nTrk==5) || (nTrk==6)){
      if((10.<massDV) && (massDV<15.)) eff =  0.01546;
      if((15.<massDV) && (massDV<20.)) eff =  0.01938;
      if((20.<massDV) && (massDV<30.)) eff =  0.0171;
      if((30.<massDV) && (massDV<50.)) eff =  0.02936;
      if((50.<massDV) && (massDV<100.)) eff = 0.04375;
   }
   if((nTrk==6) || (nTrk==7)){
      if((10.<massDV) && (massDV<15.)) eff = 0.02097;
      if((15.<massDV) && (massDV<20.)) eff = 0.03137;
      if((20.<massDV) && (massDV<30.)) eff = 0.0263;
      if((30.<massDV) && (massDV<50.)) eff = 0.04321;
      if((50.<massDV) && (massDV<100.)) eff =0.01538;

   } 
   if((7<nTrk) && (nTrk<=10)){
    if((10. <massDV) && (massDV<15.)) eff = 0.03;
    if((15. <massDV) && (massDV<20.)) eff = 0.04387;
    if((20. <massDV) && (massDV<30.)) eff = 0.05253;
    if((30. <massDV) && (massDV<50.)) eff = 0.05951;
    if((50. <massDV) && (massDV<100.)) eff = 0.05207;
    if((100.<massDV) && (massDV<200.)) eff = 0.0122;
    if((200.<massDV) && (massDV<500.)) eff = 0.0229;

 }
    if((10<nTrk) && (nTrk<=15)){
    if((10. <massDV) && (massDV<15.)) eff = 0.02161;
    if((15. <massDV) && (massDV<20.)) eff = 0.04796;
    if((20. <massDV) && (massDV<30.)) eff = 0.08033;
    if((30. <massDV) && (massDV<50.)) eff = 0.0959;
    if((50. <massDV) && (massDV<100.)) eff = 0.08688;
    if((100.<massDV) && (massDV<200.)) eff = 0.0532;
    if((200.<massDV) && (massDV<500.)) eff = 0.04351;
    if((500.<massDV) && (massDV<5000.)) eff = 0.02956;
   } 
   if((15<nTrk) && (nTrk<=20)){
    if((10. <massDV) && (massDV<15.)) eff = 0.003922;
    if((15. <massDV) && (massDV<20.)) eff = 0.02068;
    if((20. <massDV) && (massDV<30.)) eff = 0.06417;
    if((30. <massDV) && (massDV<50.)) eff = 0.114;
    if((50. <massDV) && (massDV<100.)) eff = 0.1059;
    if((100.<massDV) && (massDV<200.)) eff = 0.08409;
    if((200.<massDV) && (massDV<500.)) eff = 0.06343;
    if((500.<massDV) && (massDV<5000.)) eff = 0.04305;
 }
  if((20<nTrk) && (nTrk<=30)){
    if((20. <massDV) && (massDV<30.)) eff = 0.03704;
    if((30. <massDV) && (massDV<50.)) eff = 0.06935;
    if((50. <massDV) && (massDV<100.)) eff = 0.1007;
    if((100.<massDV) && (massDV<200.)) eff = 0.1197;
    if((200.<massDV) && (massDV<500.)) eff = 0.1048;
    if((500.<massDV) && (massDV<5000.)) eff = 0.07524;
 }
    if((30<nTrk) && (nTrk<=50)){
     if((30. <massDV) && (massDV<50.)) eff = 0.00961;
     if((50. <massDV) && (massDV<100.)) eff = 0.0753;
     if((100.<massDV) && (massDV<200.)) eff = 0.138;
     if((200.<massDV) && (massDV<500.)) eff = 0.154;
     if((500.<massDV) && (massDV<5000.)) eff = 0.134;
    }
    if((50<nTrk) && (nTrk<=200)){
     if((100.<massDV) && (massDV<200.)) eff = 0.1051;
     if((200.<massDV) && (massDV<500.)) eff = 0.2146;
     if((500.<massDV) && (massDV<5000.)) eff = 0.2051;
    }
 }
  return eff;

 }

}