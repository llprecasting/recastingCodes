/*

  CMS EXO-19-010 disappearing track search 139 fb^-1   
  -- arXiv:2004.05153  
  Recast By Mark Goodsell (goodsell@lpthe.jussieu.fr) 
  -- arXiv:2106.XXXXX

  Reuse of this code is permitted provided credit is given!
  
*/ 



#include "DT_CMS.h"
#include "include/isolation.h"

void get_which_hits(const std::vector<double> &radii, std::vector<bool> &outhits, double maxR, double Rdec, double Rprod)
{
  //std::vector<bool> outhits;
  for(auto discR : radii)
  {
    if(discR > maxR) return;

    if(Rprod > discR)
    {
      outhits.push_back(false);
      continue;
    }

    if(Rdec < discR)
    {
      outhits.push_back(false);
      continue;
    }
    
    outhits.push_back(true);

  }

  
}

std::vector<bool> get_truth_hits(const HEP::Particle* p)
{
  // This function to determine how many tracker layers are hit

  // strategy: Determine truth level and then apply per-hit efficiencies
  // Start with pixel detector and work outwards: TIB, TID, TEC, TOB
  // Pixel detector: https://lss.fnal.gov/archive/design/fermilab-design-2012-02.pdf figure 1.2, 2.1, table 2.1
  // Tracker: 1405.6569

  std::vector<bool> outhits;
  outhits.clear();

  
  double abseta=p->abseta();

  // don't bother with large eta
  if(abseta > 2.5)
  {
    //outhits.push_back(false);
    return outhits;
  }

  HEP::P4 vstart = p->prod_vertex();

  double ppz = fabs(p->mom().pz());
  double ppT=p->pT();
  
  double ootantheta=0.0;
  if(ppz > 0.0)
  {
    ootantheta = ppz/ppT; // nb for abseta > 2.5 we don't worry about zeroes here
  }

  double vdecz, vdecR;
  if(p->does_decay())
  {
    HEP::P4 vdec = p->decay_vertex();
    vdecz=fabs(vdec.pz());
    vdecR=vdec.pT();
  }
  else
  {
    vdecR=1.0e6;
    vdecz=vdecR*ootantheta;
  }

  // for now assume tracks start at z=0
  double vstartz=fabs(vstart.pz());
  double vstartR=vstart.pT();
  

  

  // Pixel
  // https://arxiv.org/pdf/0911.5434.pdf for 
  // Pixel barrel layers at e = 3, 6.8, 10.2, 16 cm
  // pixel z extends to |z| = 53.3/2 cm
  // This means |eta| < 1.3 is entirely inside the pixel barrel
  // Pre-upgrade: The endcap disks, extending from 6 to 15 cm in radius, are placed at z = ±35.5 cm and z = ±48.5 cm.
  // Post-upgrade: three endcap disks, see figure 2.1 in https://lss.fnal.gov/archive/design/fermilab-design-2012-02.pdf
  // let us guess that these are at 30, 40, 50 cm

  //1710.03842: Theradii of the four barrel layers are 2.9 cm, 6.8 cm, 10.9 cm, and 16.0 cm. 
  //The three forward disks arepositioned alongzat 3.2 cm, 3.9 cm, and 4.8 cm. 
  // Each disk is composed of two rings of moduleswith average radii of 12.8 cm and 7.8 cm.

  //static vector<double> pixelbarrelR= {30.0,68.0,102.0,160.0};
  static vector<double> pixelbarrelR= {29.0,68.0,109.0,160.0};  


  //static vector<double> pixelForwardZ = {300.0,400.0,500.0};
  static vector<double> pixelForwardZ = {320.0,390.0,480.0};


  // cannot escape inner barrel of pixel altogether
  double maxpixelR;
  if(ootantheta > 1)
  {
    maxpixelR=266.5/ootantheta;
  }
  else
  {
    maxpixelR=200.0;
  }

  get_which_hits(pixelbarrelR, outhits, maxpixelR, vdecR, vstartR);

  // check if we need pixel endcaps
  if(abseta > 1.3) 
  {
    
    double maxpixelZ=150.0*ootantheta;
    
    get_which_hits(pixelForwardZ, outhits, maxpixelZ, vdecz, vstartz);

  }
  
  if(outhits.size() < 4) 
  {
    //std::cout << "Track shorter than 4 pixel hits! " << abseta << std::endl;
    // Have I done something wrong for this case? Here I will pad the tracks with "false"
    // because we first check if there are any missing pixel hits
    for(int i=0; i < 4-outhits.size(); i++)
    {
      outhits.push_back(false);
    }

  }




  // Now for tracker barrel, see CMS 2008 report for some of these
  /*
  The Tracker Inner Barrel (TIB)
and Disks (TID) cover r < 55 cm and | z | < 118 cm, and are composed of four barrel layers,
supplemented by three disks at each end. The Tracker Outer Barrel (TOB) covers r > 55 cm and
| z | < 118 cm and consists of six barrel layers. The Tracker EndCaps (TEC) cover the region 124 <
| z | < 282 cm. Each TEC is composed of nine disks, each containing up to seven concentric
rings of silicon strip modules, yielding a range of resolutions similar to that of the TOB.
*/

  
  //static vector<double> TIBR= {230.0, 300.0, 400.0,500.0};
  // See CMS 2008 report page 64, these exetend -700 to 700 mm
  static vector<double> TIBR= {255.0, 339.0, 418.5,498.0};

  // The TID± are assemblies of three disks placed in z between ±800 mm and ±900 mm. The
  // disks are identical and each one consists of three rings which span the radius from roughly 200 mm
  // to 500 mm.
  // Here I use values inferred from Figure 3.34 in the CMS 2008 report
  static vector<double> TIDZ= {775.0,900.0,1025.0};

  // Close!!
  //static vector<double> TOBR = {600.0,700.0,800.0,890.0,980.0,1060.0};
  // See CMS 2008 report page 68
  static vector<double> TOBR = {608.0,692.0,780.0,868.0,965.0,1080.0};

  // Inferred from Figure 3.34 in the CMS 2008 report
  static vector<double> TECZ= {1250.0,1400.0,1550.0,1700.0,1950.0,2000.0,2225.0,2450.0,2700.0};

  // inner barrel
  double maxIBR;

  // inner barrel up to |z| = 580mm, r < 550
  //if(ootantheta < 1.05)
  // inner barrel up to |z| = 700mm, r < 550 -> 1/tanbeta < 700/550 = 1.27
  if(ootantheta < 1.27)
  {
    maxIBR=550.0;
    get_which_hits(TIBR, outhits, maxIBR, vdecR, vstartR);
  }
  else
  {
    maxIBR=550.0/ootantheta;
    get_which_hits(TIBR, outhits, maxIBR, vdecR, vstartR);
    //outer inner endcap, only possible if ootantheta > 1.05
    
    //double maxIBZ=580.0*ootantheta;
    double maxIBZ=700.0*ootantheta;
    get_which_hits(TIDZ, outhits, maxIBZ, vdecz, vstartz);
  }
  

  // now outer barrel, r > 55 cm and | z | < 118 cm 
  if(ootantheta < 2.145)
  {
    double maxOBR=1200.0;
    get_which_hits(TOBR, outhits, maxOBR, vdecR, vstartR);

  }
  else
  {
    double maxOBR=1180.0/ootantheta;
    get_which_hits(TOBR, outhits, maxOBR, vdecR, vstartR);

    // Now check the endcaps

    double maxTEZ=1180.0*ootantheta;
    get_which_hits(TECZ, outhits, maxTEZ, vdecz, vstartz);
  }




  // Endcaps




  /*
  std::cout << "Track with abseta = " << abseta << " and decz, decR = " << vdecz<< ", " << vdecR << std::endl;
  for(auto hit : outhits)
  {
    if(hit) 
    { 
      std::cout << 1.0 << "," ;
    }
    else
    {
      std::cout << 0.0 << "," ;
    }
  }
  std::cout << std::endl;
  */
  return outhits;
}


namespace HEP {


class charged_track  {
   private:
    
   public:
    const Particle* p; // in future change this to a vector or something else
    std::vector<bool> _hits;
    
    charged_track() { }

    charged_track(Particle* q)
    {
      p = q;
      _hits = get_truth_hits(p);
      
    };

    charged_track(Particle* q, std::mt19937 &engine, std::uniform_real_distribution<double> &rd)
    {
      p = q;
      std::vector<bool> newhits = get_truth_hits(p);

      // Now model efficiencies: first a per-hit efficiency
      // 94.5 percent is really for the earlier data; the later data claims 99% although I don't know whether to believe it
      int nhits=0;
      for(auto hit : newhits)
      {
        if(hit)
        {
          double prob = rd(engine);
          //cout << prob << ", " ;
          if(prob < 0.945)
          {
            _hits.push_back(true);
            nhits++;
            //std::cout << "1,";
            continue;
          }
          
        }
        _hits.push_back(false);
        //std::cout << "0,";
      }
      //cout << std::endl;

      // Next: we know that tracks with a smaller number of total hits 
      // see 1405.6569 page 21, "As a consequence, weaker selection criteria can be applied for tracks having many hit layers"
      // Apply 70% chance of reconstructing if nhits=4, 80% if nhits=5, 90% if 6 100% if > 6
      
      if(nhits < 7)
      {
        double recoprob=1.0;
        if(nhits == 4) recoprob=0.7;
        else if (nhits == 5) recoprob = 0.8;
        else if (nhits == 6) recoprob = 0.9;
        
        if(rd(engine) > recoprob)
        {
          _hits.clear();
          //cout << "Cleared hits! " << _hits.size() << std::endl;

        }

      }


      //std::cout << std::endl;
    }

    
    double pT() { return p->pT();}

    double E() { return p->E(); }

    double abseta() { return p->abseta(); }

    double eta() {return p->eta(); }

    double phi() { return p->phi(); }
    const HEP::P4& mom() const { return p->mom(); }
    int pid() {return p->pid(); }

    int abspid() {return p->abspid();} 

    int get_3Q() {return p->get_3Q();}

};

}


void DT_CMS::init() {
   // cout << "initialising analysis" << endl;


  cout << "--------------------------------------------------" <<std::endl;
  cout << "-- CMS disappearing track search 139 fb^-1      --" <<std::endl;
  cout << "-- arXiv:2004.05153                             --" << std::endl;
  cout << "-- By Mark Goodsell (goodsell@lpthe.jussieu.fr) --" << std::endl;
  cout << "--------------------------------------------------" <<std::endl;

  //    std::string preliminary[] ={"Njets25 >=2","1 signal lepton","Second baseline lepton veto","mT>50","ET>180","Njets <=3","2 bjets","mbb>50","ET>240","mbb"};
    //const std::vector<std::string> allSRs={"SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B"};

    const std::vector<std::string> allSRs={"SR3_2015","SR3_2016","SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B"};

    const std::vector<std::string> allSRs2={"SR3_2016","SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B"};


    const std::vector<std::string> allSR1={"SR1_2017","SR1_2018A","SR1_2018B"};
    const std::vector<std::string> allSR2={"SR2_2017","SR2_2018A","SR2_2018B"};
    const std::vector<std::string> allSR3={"SR3_2017","SR3_2018A","SR3_2018B"};

    const std::vector<std::string> all2017={"SR1_2017","SR2_2017","SR3_2017"};
    const std::vector<std::string> all2018A={"SR1_2018A","SR2_2018A","SR3_2018A"};
    const std::vector<std::string> all2018B={"SR1_2018B","SR2_2018B","SR3_2018B"};
    
    for(std::string region : allSRs)
      {
	AddRegionSelection(region);
      }



  /* 2016
total,116.0,1.0,-1.0
trigger,21.8,0.2,-0.2
passes $p_{\text{T}}^{\text{miss}}$ filters,21.8,0.2,-0.2
"$p_{\text{T}}^{\text{miss}} > 100\,\text{GeV}{}$",21.3,0.2,-0.2
"$\geq 1$ jet with $p_{\text{T}} > 110\,\text{GeV}{}$",17.3,0.2,-0.2
$\Delta\phi_{\text{max}} < 2.5$,14.9,0.2,-0.2
"$|\Delta\phi(\mbox{leading jet}, \vec{p}_{\text{T}}^{\text{miss}})| > 0.5$",14.9,0.2,-0.2
$\geq 1$ track with $|\eta| < 2.1$,14.8,0.2,-0.2
"$\geq 1$ track with $p_{\text{T}} > 55\,\text{GeV}{}$",11.5,0.2,-0.2
$\geq 1$ track passing fiducial selections,8.4,0.1,-0.1
$\geq 1$ track with $\geq 3$ pixel hits,7.2,0.1,-0.1
$\geq 1$ track with $\geq 7$ tracker hits,5.8,0.1,-0.1
$\geq 1$ track with no missing inner hits,5.8,0.1,-0.1
$\geq 1$ track with no missing middle hits,4.8,0.1,-0.1
$\geq 1$ track with relative track isolation $< 5\%$,3.8,0.1,-0.1
"$\geq 1$ track with $|d_0| < 0.02\,\text{cm}{}$",3.8,0.1,-0.1
"$\geq 1$ track with $|z_0| < 0.5\,\text{cm}{}$",3.8,0.1,-0.1
"$\geq 1$ track with $\Delta R (\text{track}, \text{jet}) > 0.5$",3.8,0.1,-0.1
"$\geq 1$ track with $\Delta R (\text{track}, \text{electron}) > 0.15$",3.53,0.09,-0.09
"$\geq 1$ track with $\Delta R (\text{track}, \text{muon}) > 0.15$",3.46,0.09,-0.09
"$\geq 1$ track with $\Delta R (\text{track}, \tau_{\text{h}}) > 0.15$",3.46,0.09,-0.09
"$\geq 1$ track with $E_{\text{calo}} < 10\,\text{GeV}{}$",3.37,0.09,-0.09
$\geq 1$ track with $\geq 3$ missing outer hits,2.01,0.07,-0.07


  */
   
    // Preselection 
  /*  2017
    total,237.6,1.1,-1.1
trigger,45.2,0.5,-0.5
passes $p_{\text{T}}^{\text{miss}}$ filters,45.2,0.5,-0.5
"$p_{\text{T}}^{\text{miss}} > 120\,\text{GeV}{}$",43.6,0.5,-0.5
"$\geq 1$ jet with $p_{\text{T}} > 110\,\text{GeV}{}$ and $|\eta| < 2.4$",32.7,0.4,-0.4
"==0 pairs of jets with $\Delta\phi_{\text{jet, jet}} > 2.5$",28.3,0.4,-0.4
"$|\Delta\phi(\mbox{leading jet}, \vec{p}_{\text{T}}^{\text{miss}})| > 0.5$",28.3,0.4,-0.4
$\geq 1$ track with $|\eta| < 2.1$,28.1,0.4,-0.4
"$\geq 1$ track with $p_{\text{T}} > 55\,\text{GeV}{}$",22.8,0.3,-0.3
$\geq 1$ track passing fiducial selections,16.5,0.3,-0.3
$\geq 1$ track with $\geq 4$ pixel hits,11.9,0.2,-0.2
$\geq 1$ track with no missing inner hits,11.8,0.2,-0.2
$\geq 1$ track with no missing middle hits,11.1,0.2,-0.2
$\geq 1$ track with relative track isolation $< 5\%$,9.03,0.21,-0.21
"$\geq 1$ track with $|d_{\text{xy}}| < 0.02\,\text{cm}{}$",9.03,0.21,-0.21
"$\geq 1$ track with $|d_z| < 0.5\,\text{cm}{}$",9.02,0.21,-0.21
"$\geq 1$ track with $\Delta R (\text{track}, \text{jet}) > 0.5$",8.87,0.21,-0.21
"$\geq 1$ track with $\Delta R (\text{track}, \text{electron}) > 0.15$",8.35,0.2,-0.2
"$\geq 1$ track with $\Delta R (\text{track}, \text{muon}) > 0.15$",8.15,0.2,-0.2
"$\geq 1$ track with $\Delta R (\text{track}, \tau_{\text{h}}) > 0.15$",8.15,0.2,-0.2
"$\geq 1$ track with $E_{\text{calo}} < 10\,\text{GeV}{}$",8.01,0.2,-0.2
$\geq 1$ track with $\geq 3$ missing outer hits,5.34,0.16,-0.16
    */

    AddCut("MET2015","SR3_2015");
    AddCut("MET",allSRs2);

    AddCut("OneGoodJet",allSRs);

    AddCut("JetAngleCuts",allSRs);
    AddCut("|Delta phi(leading jet, pTmiss)| > 0.5",allSRs);

    AddCut(">=1 track with |eta| < 2.1",allSRs);
    AddCut(">=1 track with p_T > 55",allSRs);
    AddCut(">=1 track passing fiducial selections",allSRs);
    AddCut(">=1 track with > 3 or 4 pixel hits",allSRs);
    AddCut(">=1 track with no missing inner/middle hits",allSRs);
    AddCut(">=1 track with relative isolation",allSRs);
    AddCut(">=1 track with d0 < 0.2 mm",allSRs);
    AddCut(">=1 track with dz < 5.0 mm",allSRs);
    AddCut(">=1 track with DR(track, jet) > 0.5",allSRs);
    AddCut(">=1 track with DR(track, electron) > 0.15",allSRs);
    AddCut(">=1 track with DR(track, muon) > 0.15",allSRs);
    AddCut(">=1 track with DR(track, tau) > 0.15",allSRs);
    AddCut(">=1 track with Ecalo < 10",allSRs);
    AddCut(">=1 track with >=3 missing outer hits",allSRs);

    AddCut("HEMVeto",all2018B);
    
    //AddCut("BadAngles2017",all2017);
    //AddCut("BadAngles2018",all2018A);  
    //AddCut("BadAngles2018",all2018B);  

    //AddCut("4 Layers",allSR1);
    //AddCut("5 Layers",allSR2);
    //AddCut(">5 Layers",allSR3);
    
    AddCut(">5 Layers 2015","SR3_2015");
    AddCut(">5 Layers 2016","SR3_2016");

    AddCut("4 Layers 2017","SR1_2017");
    AddCut("4 Layers 2018A","SR1_2018A");
    AddCut("4 Layers 2018B","SR1_2018B");

    AddCut("5 Layers 2017","SR2_2017");
    AddCut("5 Layers 2018A","SR2_2018A");
    AddCut("5 Layers 2018B","SR2_2018B");

    AddCut(">5 Layers 2017","SR3_2017");
    AddCut(">5 Layers 2018A","SR3_2018A");
    AddCut(">5 Layers 2018B","SR3_2018B");
  
    vector<double> METedges={10,30,50,70,90,110,130,150,200,250,300,350,400,1000.0};
    
    AddYodaHisto1D("MET",METedges);
    AddYodaHisto1D("pTJ1",100,10.0,2000.0);
    AddYodaHisto1D("pTTrack1",100,10.0,2000.0);
    vector<double> DRedges={0.0,0.5,1.0,100.0};
    AddYodaHisto1D("CharginoDeltaR",DRedges);
   // AddYodaHisto1D("pTJ2",myedges);
//cout << "initialised analysis" << endl;
};

DT_CMS::~DT_CMS() { 
  //delete Event;

};


void DT_CMS::Execute(std::mt19937 &engine) {
  
  //std::vector<HEP::Particle*> Electrons, Muons, BaselineElectrons, BaselineMuons, BaselineLeptons, SoftElectrons, SoftMuons;

  //std::vector<HEP::Jet*> SoftJets, SignalJets, bJets, Jets25;

  // Now do the MET calculation from *visible* objects
  HEP::P4 pTmiss,pTmissNoMu;

  std::vector<HEP::Particle*> charginos;
// Count signal electrons, muons, jets etc 
//cout << "Event has " << Event.electrons().size() << "electrons and " << Event.jets().size() << " jets " << endl;
//for(auto electron : Event->electrons())

static std::uniform_real_distribution<double> rd(0.0,1.0);
bool chargino50=false;

std::vector<HEP::charged_track*> tracks;

HEP::P4 muonp;

  for(auto muon : Event->muons())
  {
    muonp+=muon->mom();

  }

double maxtrackpt=0.0;
double chpt=0.0;

//std::vector<HEP::Particle*> new_fake_muons;
  for(auto chargino : Event->HSCPs())
  {
    //double d0Calc(const HEP::P4 &xv, const HEP::P4 &p);
    //double absd0 = d0Calc(chargino->prod_vertex(),chargino->mom());
    //double absdz = dzCalc(chargino->prod_vertex(),chargino->mom());
    //if((chargino->abseta() < 2.5) && (absd0 < 0.2) && (absdz < 5.0))  // distances in mm
    //cout << "Got a chargino: " << std::endl;
    if((chargino->decay_vertex().pT() > 7.0e3) || (fabs(chargino->decay_vertex().pz()) > 1.1e4))
    {
      // gets reconstructed as a muon
      muonp+=chargino->mom();
      //new_fake_muons.push_back(chargino);
    } 
    chpt=chargino->pT();
    if(chpt > maxtrackpt) maxtrackpt=chpt;
    
    if(chargino->abseta() < 2.5)
            {
              HEP::charged_track* newtrack = new HEP::charged_track(chargino,engine,rd);
              tracks.push_back(newtrack);
              //  charginos.push_back(chargino);
              if((!chargino50) && (chargino->pT() >50.0)) chargino50=true;
              

            }

  }

  //add charged hadrons to our "charginos"; don't bother with leptons as pileup doesn't produce many
  
  for(auto chargino : Event->get_charged_hadrons())
  {
    if((chargino->decay_vertex().pT() > 7.0e3) || (fabs(chargino->decay_vertex().pz()) > 1.1e4))
    {
      // gets reconstructed as a muon
      muonp+=chargino->mom();
     // new_fake_muons.push_back(chargino);
    } 

    if(chargino->abseta() < 2.5)
            {
              HEP::charged_track* newtrack = new HEP::charged_track(chargino,engine,rd);
              tracks.push_back(newtrack);
              //  charginos.push_back(chargino);
                if((!chargino50) && (chargino->pT() >50.0)) chargino50=true;
            }
  }
  
  

  

  pTmiss=Event->missingmom();
  pTmissNoMu=pTmiss+muonp;
  
  double MET=Event->met();
  double METnomu=pTmissNoMu.pT();


  FillYodaHisto1D("MET",MET);
  FillYodaHisto1D("pTTrack1",maxtrackpt);
  if( Event->jets().size() > 0)
  {
    FillYodaHisto1D("pTJ1", Event->jets()[0]->pT());

  }

  bool passMET=false;
  bool passMET2015=false;

  if((chargino50) && (MET > 105.0)) passMET=true;
  if(MET > 120.0) passMET = true;
  if(METnomu < 120.0) 
  {
    passMET = false;
  }
  else
  {
    passMET = true;
  }

  if((chargino50) && (MET > 75.0) ) passMET2015=true;
  if((MET > 90.0) || (METnomu > 90.0)) passMET2015 = true;
  if(MET < 100.0) passMET2015 = false;

  
  double eff=1.0;
  //  1903.06078 figure 5, since the first cut on MET is the trigger. It uses HLT with 90 GeV or 120 GeV
  if(passMET2015)
  {

    double tMET=MET;
    if (METnomu > MET) tMET=METnomu;
    if(tMET < 90.0) 
    {
      eff=0.0;
    } 
    else if (tMET < 250.0)
    {
      eff = (tMET-50.0)/200.0;
    }
    else
    {
      eff=1.0;
    }
    
    if(rd(engine) > eff) passMET2015=false;
  }

  if(passMET)
  {
    // Use METnomu?

    double tMET=MET;
    if (METnomu > MET) tMET=METnomu;
    
    if(tMET < 120.0) 
    {
      eff=0.0;
    } 
    else if (tMET < 280.0)
    {
      eff = (tMET-120.0)/160.0;
    }
    else
    {
      eff=1.0;
    }

    if(rd(engine) > eff) passMET=false;
  }
 


  ApplyCut(passMET,"MET");
  ApplyCut(passMET2015,"MET2015");

  // need to clean up all the tracks I've new'd
  if((!passMET) || (!passMET2015)) {for(auto track: tracks) {delete track;}; return;} ;
  bool onegoodjet = false;
  bool JetMetDeltaphi=true;
  bool passJetAngleCuts = true;


  std::vector<HEP::Jet*> Jets = Event->jets();
  std::vector<HEP::Jet*> GoodJets;
  std::vector<HEP::Particle*> GoodCharginos;

 
  if(Jets.size() ==0 ) {for(auto track: tracks) {delete track;}; return;} ;
  if(Jets[0]->mom().deltaPhi(pTmiss) < 0.5) JetMetDeltaphi = false;

  filterPhaseSpace(Jets,30.0,2.4); 

  for(int m = 0; m< Jets.size() ; m++)
  {
    HEP::Jet* tjet=Jets[m];
   // if(tjet->pT() < 20.0) continue;
    if((!onegoodjet) && (tjet->abseta()< 2.4) && (tjet->pT() > 110.0)) 
    {
      onegoodjet=true;
      //break;
    }
    if(tjet->pT() > 30.0) 
    { 
      GoodJets.push_back(tjet);

    }
    for(int n=m+1; n < Jets.size(); n++)
    {
      //if(Jets[n]->pT() < 20.0) continue;
      if(fabs(tjet->mom().deltaPhi(Jets[n]->mom())) > 2.5)
      {
        passJetAngleCuts = false;
        break;
      }

    }

  }

  

  ApplyCut(onegoodjet,"OneGoodJet");
  if(!onegoodjet)  {for(auto track: tracks) {delete track;}; return;} ;
  

    ApplyCut(passJetAngleCuts,"JetAngleCuts");
  ApplyCut(JetMetDeltaphi,"|Delta phi(leading jet, pTmiss)| > 0.5");
  if((!passJetAngleCuts) || (!JetMetDeltaphi))  {for(auto track: tracks) {delete track;}; return;} ;

      


// ought to do this for taus as well ...

if(tracks.size() < 1)  return;
  bool good2015_3=false;
  bool good2016_3=false;
  bool good2017_1=false;
  bool good2017_2=false;
  bool good2017_3=false;
  bool good2018A_1=false;
  bool good2018A_2=false;
  bool good2018A_3=false;
  bool good2018B_1=false;
  bool good2018B_2=false;
  bool good2018B_3=false;

  bool tracketa=false;
  bool trackpT=false;
  bool trackfiducial=false;
  bool trackpixel=false;
  bool trackinnermiddle=false; 
  bool trackisolation=false;
  bool trackd0=false;
  bool trackdz=false;
  bool trackDRjet=false;
  

  std::vector<HEP::charged_track*> tracks2;

  for(auto chargino: tracks)
  {
    if((chargino->abseta() > 2.1) ) {delete chargino; continue;};

    tracketa=true;

    if(chargino->pT() < 55.0) {delete chargino; continue;};
    trackpT=true;

    double cheta=chargino->abseta();

    // apply conditions on the track to avoid low efficiency parts of muon chamber/ECAL
    // muon chamber
    //if((decayrad > 4020.0) && (decayrad < 7380.0))
    //{
      if((cheta >0.15) && (cheta < 0.35) ) {delete chargino; continue;};
     if((cheta >1.55) && (cheta < 1.85) ) {delete chargino; continue;};
    //}
    // ECAL
    //if((decayrad > 1290.0) && (decayrad < 1790.0))
    //{
      if((cheta >1.42) && (cheta < 1.65) ) {delete chargino; continue;};

    trackfiducial=true;


    // Now we have >=1 track with > 3 or 4 pixel hits, actually this means we have to have all 4 closest hits

    if(chargino->_hits.size() < 4) {delete chargino; continue;};

    bool charginotrackpixel=true;
    
    for(int i=0; i<4; i++)
    { 
      if(!chargino->_hits[i])
      {
        charginotrackpixel=false;
        break;
      }
    }
    if(!charginotrackpixel) {delete chargino; continue;};
    trackpixel=true;
    bool charginoinnermiddle=true;
    // now check for missing inner/middle hits
    if(chargino->_hits.size() > 4)
    {
      bool foundoutside=false;
       for(int i=chargino->_hits.size() -1; i>3; i--)
       {
         if(foundoutside)
         {
           if(!chargino->_hits[i])
           {
             charginoinnermiddle=false;
             break;
           }

         }
         else
         {
           if(chargino->_hits[i])
           {
             foundoutside=true;
           }
         }
       } 

    }

    if(!charginoinnermiddle) {delete chargino; continue;};

    trackinnermiddle=true;

    double isopt = sumpTisolation(chargino->p->mom(),Event,0.3); 
    //std::cout << "pT isolation: " << isopt << " pT: " << chargino->pT() << " ratio: " << (isopt/chargino->pT())<< std::endl;
    if( (isopt/chargino->pT()) > 0.05) {delete chargino; continue;};

    trackisolation=true;
    
    double absd0 = d0Calc(chargino->p->prod_vertex(),chargino->p->mom());
    double absdz = dzCalc(chargino->p->prod_vertex(),chargino->p->mom());

    if(absd0 > 0.2) {delete chargino; continue;};
    trackd0=true;
    if(absdz > 5.0) {delete chargino; continue;};
    trackdz=true;

    tracks2.push_back(chargino);
  }

  ApplyCut(tracketa,">=1 track with |eta| < 2.1");
  ApplyCut(trackpT,">=1 track with p_T > 55");
  ApplyCut(trackfiducial,">=1 track passing fiducial selections");
  ApplyCut(trackpixel,">=1 track with > 3 or 4 pixel hits");
  ApplyCut(trackinnermiddle,">=1 track with no missing inner/middle hits");
  ApplyCut(trackisolation,">=1 track with relative isolation");
  ApplyCut(trackd0,">=1 track with d0 < 0.2 mm");
  ApplyCut(trackdz,">=1 track with dz < 5.0 mm");



if(tracks2.size() < 1) {for(auto track: tracks2) {delete track;}; return;} ;

if( Event->HSCPs().size() == 2)
   {
     FillYodaHisto1D("CharginoDeltaR",Event->HSCPs()[0]->mom().deltaR_eta(Event->HSCPs()[1]->mom()));

   }

tracks=FullRemoval(tracks2,GoodJets,0.5);
ApplyCut((tracks.size() > 0),">=1 track with DR(track, jet) > 0.5");
if(tracks.size() ==0) return;
tracks=FullRemoval(tracks,Event->electrons(),0.15);
ApplyCut((tracks.size() > 0),">=1 track with DR(track, electron) > 0.15");
if(tracks.size() ==0) return;

// This cut should remove long-lived charginos which get reconstructed as a muon
tracks=FullRemoval(tracks,Event->muons(),0.15);

std::vector<HEP::charged_track*> tracks3;
std::vector<HEP::charged_track*> fake_muons;
  for(auto chargino: tracks)
  {
    HEP::P4 vdec = chargino->p->decay_vertex();
    //if((vdec.pT() > 4000.0) || (fabs(vdec.pz()) > 6000.0 )) {delete chargino; continue;}
    //if((vdec.pT() > 7000.0) || (fabs(vdec.pz()) > 11000.0 )) {delete chargino; continue;}
    // if it decays outside muon chamber
    if((vdec.pT() > 7000.0) || (fabs(vdec.pz()) > 11000.0 ))
    {
      fake_muons.push_back(chargino);
    }
    else
    {
      tracks3.push_back(chargino);
    }
  }
  //tracks=tracks3;
  tracks=FullRemoval(tracks3,fake_muons,0.15);
  for(auto fmu : fake_muons)
  {
    delete fmu;
  }

ApplyCut((tracks.size() > 0),">=1 track with DR(track, muon) > 0.15");
if(tracks.size() ==0) return;
tracks=FullRemoval(tracks,Event->taus(),0.15);
ApplyCut((tracks.size() > 0),">=1 track with DR(track, tau) > 0.15");

if(tracks.size() ==0) return;

bool trackEcalo=false;
bool trackmissingouter=false;

  // this is according to Hart's thesis the way of doing pileup subtraction
  // we take the "Median" calorimeter energy (whatever that is) and multiply by the cone area
  

  for(auto chargino: tracks)
  {
   // if((chargino->pT() < 55.0) || (chargino->abseta() > 2.1) ) continue;

    //double isopt = sumPTisolation(chargino->mom(),Event,0.3)+chargino.pT(); // since in sumPTisolation we subtract the track pT, but the charginos aren't included!

    // this is according to Hart's thesis the way of doing pileup subtraction
    // we take the "Median" calorimeter energy and multiply by the cone area
    /*
    double rhofactor=(Event->Average_rho)*coneAreaApprox(chargino->p->mom(),0.5);
    double isocalo = sumCaloE(chargino->p->mom(),Event,0.5)-rhofactor;
    
    std::cout << "isoCaloE: " << isocalo << " rho: " << rhofactor << std::endl;

    if(isocalo > 10.0) continue;
    */

    //double rhofactor=0.0;
    //double isocalo = sumCaloE(chargino->p->mom(),Event,0.5)-rhofactor;
    //if (isocalo > 10.0) continue;
    trackEcalo=true;
    

    
    

 
      double chphi=chargino->phi();
     double cheta=chargino->eta();
    double acheta=fabs(cheta);

    

    bool chtrackmissingouter=true;
    int nhits=(chargino->_hits.size())-1;
    for(int i=0; i<3; i++) 
       {
         if(chargino->_hits[nhits-i])  // 3 outermost hits must be missing
         {
           chtrackmissingouter=false;
           break;
         }

       }
    if(!chtrackmissingouter) continue;
    trackmissingouter=true;
    bool good2017 = true;
    bool good2018 = true;
    
    

    // apply conditions on the track to avoid low efficiency parts of muon chamber/ECAL
    // muon chamber
    //if((decayrad > 4020.0) && (decayrad < 7380.0))
    //{
   //   if((cheta >0.15) && (cheta < 0.35) ) continue;
    // if((cheta >1.55) && (cheta < 1.85) ) continue;
    //}
    // ECAL
    //if((decayrad > 1290.0) && (decayrad < 1790.0))
    //{
     // if((cheta >1.42) && (cheta < 1.65) ) continue;
    //}



    if((cheta >0.0) && (cheta < 1.42) && (chphi < 3.142) && (chphi > 2.7) ) good2017=false;
    if((cheta >0.0) && (cheta < 1.42) && (chphi < 0.8) && (chphi > 0.4) ) good2018=false;

    // Now we've already selected tracks that have no missing hits from start to end, so just count layers

    int nlayers=0;

    for(auto hit : chargino->_hits)
    {
      if(hit) nlayers++;
    }

    if(nlayers == 4)
    {
      if(good2017) good2017_1=true;
      if(good2018)
        {
          good2018A_1=true;
          good2018B_1=true;
        }
    }
    else if(nlayers == 5)
    {
      if(good2017) good2017_2=true;
        if(good2018)
        {
          good2018A_2=true;
          good2018B_2=true;
        }

    }
    else // nlayers > 5
    {
        if(nlayers > 7)
        {
          good2015_3=true;
          good2016_3=true;
        }
        if(good2017) good2017_3=true;
       
        if(good2018)
        {
          good2018A_3=true;
          good2018B_3=true;
        }
    }

  }

ApplyCut(trackEcalo,">=1 track with Ecalo < 10");
ApplyCut(trackmissingouter,">=1 track with >=3 missing outer hits"); 


double phipTmiss = pTmiss.phi(); 
bool GoodpTPhi=true;
if((phipTmiss > -1.6) && (phipTmiss <-0.6)) GoodpTPhi=false;

ApplyCut(GoodpTPhi,"HEMVeto");
//if(!GoodpTPhi) return;



ApplyCut(good2017_1,"4 Layers 2017");
ApplyCut(good2018A_1,"4 Layers 2018A");
ApplyCut(good2018B_1,"4 Layers 2018B");

ApplyCut(good2017_2,"5 Layers 2017");
ApplyCut(good2018A_2,"5 Layers 2018A");
ApplyCut(good2018B_2,"5 Layers 2018B");

ApplyCut(good2015_3,">5 Layers 2015");
ApplyCut(good2016_3,">5 Layers 2016");
ApplyCut(good2017_3,">5 Layers 2017");
ApplyCut(good2018A_3,">5 Layers 2018A");
ApplyCut(good2018B_3,">5 Layers 2018B");
for(auto track: tracks) {delete track;};
   
    

};
