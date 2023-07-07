// Pythia8 code to test CMS analysis for emerging jets EJs (1810.10069), giving pythia_card with details of model.
// To reproduce CMS results, use cardEJcms.cmnd
// This implemets all selection sets in Table 3, aplying cuts at track- and event-level. 

#include <Pythia8/Pythia.h> 
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <algorithm>

using namespace Pythia8;

bool hadronization_ok (Event event) {
  for (int i = 0; i < event.size(); ++i){
    if(event[i].isFinal() && event[i].id()==4900101){
      return false;
    }
  }
  return true; 
}

double cmedia (vector<double> &v){
  int size=v.size();
  sort(v.begin(), v.end());
  if(size % 2 != 0){
    return v[size/2.];
  } 
  else{
    return (v[(size-1)/2.] + v[size/2.])/2.0;
  }
}
///Compute transverse impact parameter of track
double d0t ( Particle* part){
  return abs(part->xProd()*part->py() - part->yProd()*part->px())/part->pT();
}
///Compute lomgitudinal impact parameter
double d0z ( Particle* part){
  Vec4 a = part->vProd();
  Vec4 b = part->p();
  double dt = abs(part->xProd()*part->py() - part->yProd()*part->px())/part->pT();
  double d0 = sqrt(dot3(cross3(a,b),cross3(a,b))) / part->pAbs();
  double dz = sqrt(pow(d0,2)-pow(dt,2));
  return abs(dz);
}

/// Generate random number from normal distribution
double randN (double const mu, double const sigma, int const seed){  
  //std::random_device rd;
  //std::mt19937 urn(rd());
  std::mt19937 urn(seed);
  std::normal_distribution<double> norm{mu,sigma};
  double value = norm(urn);
  double posval = abs(value);
  return posval;
}
/// Generate random number from binomial distribution
double randBool (double const pv, int const seed){
  //std::random_device rd;  
  //std::mt19937 urn(rd());
  std::mt19937 urn(seed);  
  std::bernoulli_distribution flip{pv};
  bool pass = flip(urn);
  return pass;
}
/// Read track reconstruction effiency as function of radial distance from IP, from iteration 5 (taken from Fig.12 in 1405.6569)
double readEff (double const rad){
  int bins=30;
  double wbin=2;
  int nb=-1;
  for (int i=0; i<bins; i++){
    double low = i*wbin;
    double upp = low + wbin;
    if (rad >= low && rad < upp){
      nb = i;
      break;
    }
  }
  ifstream fin;
  fin.open("Iter_5.csv", ios::in);
  vector<string> row;
  string line, val;
  double row2=0;
  int nrow=0;
  while (getline (fin, line)) {
    row.clear();
    stringstream s(line);
    while (getline(s, val, ',')) {
      row.push_back(val);
    }
    row2 = stod(row[1]);
    //if (rad >= (round(row1)-1) && rad <(round(row1)+1)) break;
    if (nrow == nb){
      goto end_fun;
    }
    nrow++;
  }
  if (fin.eof( ) && nb==-1) row2=0; 
  end_fun:
  fin.close();
  return row2;
}
/// Read resolution in d0_T of track as function of momentum and pseudorapidity (taken from Fig.15 in 1405.6569)
double resolIPpt (double const xp, double const eta){
  ifstream fin;
  fin.open("fig14_d0_resolution_pT.csv", ios::in);
  vector<string> rowL, rowU;
  string line, val;
  double x1=0,x2=1,y1=0.5,y2=0.5;
  double yp=0;
  int nrow=0;
  int ptl=0, resl=0;
  if (eta>=0 && eta<0.9){
    ptl=0;
    resl=1;
  }else if(eta>=0.9 && eta<1.4){
    ptl=4;
    resl=5;
  }else if(eta>=1.4 && eta<2.5){
    ptl=8;
    resl=9;
  }
  getline(fin, line);
  while (getline (fin, line)) {
    if (nrow>=2){
      rowL.clear();
      rowL = rowU;
      rowU.clear();
    }
    stringstream ss(line);
    while (getline(ss, val, ',')) {
      if (nrow==0){
        rowL.push_back(val);
      }else if (nrow==1){
        rowU.push_back(val);
      }else if (nrow>=2){
        rowU.push_back(val);
      }
    }
    if(nrow>=1){
      x1 = stod(rowL[ptl]);
      y1 = stod(rowL[resl]);
      x2 = stod(rowU[ptl]);
      y2 = stod(rowU[resl]);
      if (xp>=x1 && xp<x2){
        yp = y1 + ((y2-y1)/(x2-x1)) * (xp - x1);
        goto end_fun;
      }
    }   
    nrow++;
  }
  if (fin.eof( ) && xp>=x2 ) yp = y1 + ((y2-y1)/(x2-x1)) * (xp - x1);
  end_fun:
  fin.close();
  return yp;
}

int main(int argc, char* argv[]) {

  Pythia pythia;
  Event& event = pythia.event;

  if (argc != 2) {
    cout << "Give input card. Usage ./main1002.x inputcard" << endl;
    return 0;
  }

  //pythia.readFile("cardEJcms.cmnd");
  pythia.readFile(argv[1]);

  
  SlowJet slowJet( -1, 0.4, 20., 2.0, 2, 1); //power, R, ptjetmin, etamax, which particles, mass

  pythia.init();

  Event trimmedevent;
  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  //int nEvent = 1000;
  int nAbort = 30;
  
  //ofstream myfile1;
  //myfile1.open("effisR_test.dat", std::ios_base::app);
  //ofstream myfile2;
  //myfile2.open("effisIt_test.dat", std::ios_base::app);

  ///cuts on variables for EJs tagging (table2) //EMJ-G(1-6) [cz0,cdn,cip,calf] 
	vector<vector<double>> EJg{{25.0,4.,0.5,0.25},{40.0,4,1.0,0.25},{40.0,20,2.5,0.25},{25.0,4,1.0,0.25},{25.0,20.0,0.5,0.25},{25.0,10.0,0.5,0.25}}; 
	///event selection sets (table3) //setNum(1-7) [HT, ptj(1-4), met, #ej, EJg(#)]
  vector<vector<double>> setG{{900,225,100,100,100,0,2,1},{900,225,100,100,100,0,2,2},{900,225,100,100,100,200,1,3},{1100,275,250,150,150,0,2,1},{1000,250,150,100,100,0,2,4},{1000,250,150,100,100,0,2,5},{1200,300,250,200,150,0,2,6}};

  bool smear=true;
  int jsize;
  double nEvs=0;
  int met2;
  double nEvs2=0;
  vector<double> Evtot(7,0);
  vector<double> Evtot2(7,0);

  int mN=0;
  double ctN=0;

  /// Begin event loop.
  int iAbort = 0;
  //int iEvent = 0;
  //while (iEvent < nEvent) {
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
    
    /*if (iEvent <3) {
      pythia.info.list();
      pythia.event.list();
    }*/
  
    slowJet.analyze( pythia.event );	

    int ndq=0;
    vector<int> dplist;
    vector<int> dtin;
    vector<double> dtj;
    vector<int> idq;
    int ndt=0;
    vector<int> idpi;
    int ndpi=0;
    int nsq=0;
    vector<int> isq;
    Vec4 vmet2 = Vec4(0.0,0.0,0.0,0.0);
    Vec4 vudm = Vec4(0.0,0.0,0.0,0.0);
    int npist=0;

    for (int i = 0; i < pythia.event.size(); ++i){

      /// Get information for dark pions.
      if((abs(event[i].id()) ==  4900111 || abs(event[i].id()) ==  4900211) && (abs(event[event[i].daughter1()].id())!=4900111 || abs(event[event[i].daughter1()].id())!=4900211)){       
        ndpi++;
        mN = event[i].m0();
        ctN = event[i].tau0();
        if (event[i].isFinal()){
          npist++;
          vudm += event[i].p();
        }
      }
      //// ##Save all final tracks from dquarks (4900101) 
      if((abs(event[i].id()) ==  4900101 && abs(event[event[i].daughter1()].id())!=4900101)){// || (abs(event[i].id()) ==  4900021 && abs(event[event[i].mother1()].id()) == 4900001)){
        Particle part = pythia.event[i];
        ndq++;
        idq.push_back(i);
        //save all final tracks coming from the LLP
        dplist = part.daughterListRecursive();
        sort(dplist.begin(), dplist.end());
        dplist.erase(unique( dplist.begin(), dplist.end() ), dplist.end());          
        for (int q=0; q < int(dplist.size()); q++) {
          if(event[dplist[q]].isFinal() && event[dplist[q]].pT()>1.0 && event[dplist[q]].eta()<2.4){// &&event[dplist[q]].isCharged() 
            dtin.push_back(dplist[q]);                              
            dtj.push_back(slowJet.jetAssignment(dplist[q]));
            ndt++;
          }  
        }
      }
      // #Save pair of dquarks
      if((abs(event[i].id()) == 1 && abs(event[event[i].mother1()].id())==4900001)){
        nsq++;
        isq.push_back(i);
      }
      // #Calculate MET
      if(event[i].isFinal() && (abs(event[i].id())==11 || abs(event[i].id())==13 || abs(event[i].id())==22)){
        if(event[i].pT()>10.0 && abs(event[i].eta())<2.5) vmet2 += event[i].p();
      }
      if(event[i].isFinal() && event[i].isHadron()  && event[i].pT()>0.5 && abs(event[i].eta())<4.8) vmet2 += event[i].p();
    }/////end event.size 
    
    //### Analyse jets
    jsize = int(slowJet.sizeJet()); 

    vector<double> mipj1(jsize,-1);
    vector<vector<double>> alpha1(jsize, vector<double> (6,-1));
    double ptjja=0;
    vector<int> ejn(6,0); 
    int jpass=0;
    vector<vector<bool>> ejp(jsize,vector<bool>(6,false));
    bool evl4j=false;
    vector<int> jntr1(jsize,0);
    vector<bool> jinb(jsize,false);
    vector<bool> evpf(7,false);
    bool ptTrg=false;

    vector<double> mipj2(jsize,-1);
    vector<vector<double>> alpha2(jsize, vector<double> (6,-1));
    double ptjja2=0;
    vector<int> ejn2(6,0);
    int jpass2=0;
    vector<vector<bool>> ejp2(jsize,vector<bool>(6,false));
    bool evl4j2=false;
    vector<int> jntr2(jsize,0);
    vector<bool> jinb2(jsize,false);
    vector<bool> evpf2(7,false);
    bool ptTrg2=false;

    vector<double> ptjet(jsize,0);
    int epass=0;
    vector<bool> dqinj(jsize,false);
    vector<bool> mipdq(jsize,false);
    vector<bool> sqinj(jsize,false);
    vector<bool> mipdq2(jsize,false);



    
    //### loop over jets
    for (int j=0; j < jsize; j++){
      vector<int> jpart;
      double ptall=0;
      double ptall2=0;
      vector<double> ptpv1(6,0);
      vector<double> ptpv2(6,0);
      vector<double> ipj;
      vector<double> ipj2;
      double esum=0;
      bool emp=false;
      bool jin=false;
      bool jin2=false;
      bool tpass;
      double z0cut,DNcut,IPcut,alphcut;

      ///### Check if dQ or q_sm is within a jet cone
      for (int m=0; m < int(idq.size()); m++){
        double dR = RRapPhi(slowJet.p(j),event[idq[m]].p());
        if (dR < 0.4) dqinj[j]=true;
      }
      for (int s=0; s < int(isq.size()); s++){
        double dR = RRapPhi(slowJet.p(j),event[isq[s]].p());
        if (dR < 0.4) sqinj[j]=true;
      }
      
      //### loop oves jet constituents
      jpart = slowJet.constituents(j);
      for (int t=0; t< int(jpart.size()); t++){
        int seed = iEvent+j+t;  
        double rprod = abs(sqrt(pow(event[jpart[t]].xProd(),2) + pow(event[jpart[t]].yProd(),2)));
        double rpcm = rprod*0.1;
        double tpt = event[jpart[t]].pT();
        double teta = abs(event[jpart[t]].eta());
        tpass=false;
        double eftrk, ipa, z0a;
        /// selecting tracks (with two approaches)
        if (event[jpart[t]].isCharged() && event[jpart[t]].isFinal() && event[jpart[t]].pT()>1.0 && abs(event[jpart[t]].eta())<2.5){ // && rprod<600
          if (rprod<600){
            eftrk = readEff(rpcm); //as rpcm in cm
            tpass = randBool(eftrk,seed);
          }
          double z0i=d0z(&event[jpart[t]]);
          double ipi=d0t(&event[jpart[t]]);
          double sigd0 = resolIPpt(tpt, teta) *0.001; //as sigma_d0 in [mum]
          //double sigd0 =0.05;
          if (smear){            
            ipa= randN(ipi+sigd0,sigd0,seed);
            z0a= randN(z0i,0.1,seed);
          }else{
            ipa = ipi;
            z0a = z0i;             
          }
          double DNa = sqrt(pow((z0a/0.1),2) + pow((ipa/sigd0),2));
          if (rprod< 102.0){
            ipj.push_back(ipa);
            for (int u=0; u<6; u++){
              z0cut = EJg[u][0];
              DNcut = EJg[u][1];
              if (z0a<z0cut && DNa < DNcut) ptpv1[u] += event[jpart[t]].pT();
            } 
            ptall += event[jpart[t]].pT();
            jntr1[j]++;
          }
          if (tpass){
            ipj2.push_back(ipa);
            for (int u=0; u<6; u++){
              z0cut = EJg[u][0];
              DNcut = EJg[u][1];
              if (z0a<z0cut && DNa < DNcut) ptpv2[u] += event[jpart[t]].pT();
            }
            ptall2 += event[jpart[t]].pT();
            jntr2[j]++;
          }  
        }
        if(abs(event[jpart[t]].id())==11 || event[jpart[t]].id()==22) esum += event[jpart[t]].e();
      }
      ///# claculate media(IP_trk) and alpha variables
      if (ptall>1.){
        mipj1[j] = cmedia(ipj);
        for (int g=0; g<6; g++){
        	alpha1[j][g] = ptpv1[g] / ptall;
        }
        jin = true;        
      }           
      if (ptall2>1.){
        mipj2[j] = cmedia(ipj2);
        for (int g=0; g<6; g++){
        	alpha2[j][g] = ptpv2[g] / ptall2;
        }
        jin2 = true;       
      }      

      ptjet[j] = slowJet.pT(j);

      // check energy of jets from electrons and photons is less 90%
      if (esum < (0.9*slowJet.p(j).e())){
        epass++;
        emp=true;
      }
      // selecting jets and EJs
      if (jin && emp){
        jpass++;
        jinb[j]=true;
        ptjja += ptjet[j];
        // look if jet is EJ, passing IP and alpha requirements
        for (int f=0; f<6; f++){
        	IPcut = EJg[f][2];
        	alphcut = EJg[f][3];
        	if (mipj1[j] > IPcut && alpha1[j][f] < alphcut){
          	ejn[f]++;
          	ejp[j][f]=true;
        	}
        }
      } 
      if (jin2 && emp){
        jpass2++;
        jinb2[j]=true;
        ptjja2 += ptjet[j];
        for (int f=0; f<6; f++){
        	IPcut = EJg[f][2];
        	alphcut = EJg[f][3];
        	if (mipj2[j] > IPcut && alpha2[j][f] < alphcut){
          	ejn2[f]++;
          	ejp2[j][f]=true;
        	}
        }
      }
    }//## end loop on jets

////////////

    vmet2 = vmet2 - vudm;
    met2 = vmet2.pT();

    if (ptjja > 900.) ptTrg=true;
    if (ptjja2 > 900.) ptTrg2=true;

    /// #sort jets according ptjet
    vector<int> V(jsize);
    std::iota(V.begin(),V.end(),0); 
    sort( V.begin(),V.end(), [&](int i,int j){return ptjet[i]>ptjet[j];} );

    vector<int> Vi, Vi2;
    for (int k=0; k<jsize; k++){
      if (jinb[V[k]]) Vi.push_back(V[k]);
      if (jinb2[V[k]]) Vi2.push_back(V[k]); 
    }
    
    ///## Select first 4 hardest jets and apply event requirements
    if (jpass>=4 && !evl4j){
      evl4j=true;
      for (int h=0; h<7; h++){
        int ejnf=0;
        int ptj4=0;
        int gej=setG[h][7]-1;
        int cej=setG[h][6];
        for (int w=0; w<4; w++){
          if (ejp[Vi[w]][gej]) ejnf++;
          ptj4 += ptjet[Vi[w]];
        }
        if (ptj4>setG[h][0]  && ptjet[Vi[0]]>setG[h][1] && ptjet[Vi[1]]>setG[h][2] && ptjet[Vi[2]]>setG[h][3] && ptjet[Vi[3]]>setG[h][4]){
          if (ejnf >=cej && met2>=setG[h][5]) evpf[h]=true;
        }          
      }
    }
    if (jpass2>=4 && !evl4j2){
      evl4j2=true;
      for (int e=0; e<7; e++){
        int ejnf2=0;
        int ptj24=0;
        int gej2=setG[e][7]-1;
        int cej2=setG[e][6];
        for (int w=0; w<4; w++){
          if (ejp2[Vi2[w]][gej2]) ejnf2++;
          ptj24 += ptjet[Vi2[w]];
        }
        if (ptj24 >setG[e][0]  && ptjet[Vi2[0]]>setG[e][1] && ptjet[Vi2[1]]>setG[e][2] && ptjet[Vi2[2]]>setG[e][3] && ptjet[Vi2[3]]>setG[e][4]){
          if (ejnf2 >=cej2 && met2 >=setG[e][5]) evpf2[e]=true;
        }          
      }
    }  
/////////////
    if (ptTrg && evl4j) nEvs++;
    if (ptTrg2 && evl4j2) nEvs2++;

    for (int f=0; f<7; f++){
      if (evpf[f]) Evtot[f]++;
      if (evpf2[f]) Evtot2[f]++;
    }

  }///end event loop

  cout.precision(4);
  cout << "total efficiency EJs, sn7 (r<102mm)" << Evtot[6]/nEvs <<  endl;
  cout << "total efficiency EJs, sn7 (it5)" << Evtot2[6]/nEvs2 <<  endl;
  cout << "m_dp: "<<mN<<" ct_dp: "<< int(ctN) << endl;


  //myfile1 << mS<<" "<<mN <<" "<< ctN<<" "<< nEvs <<" "<< Evtot[0] <<" "<< Evtot[1] <<" "<< Evtot[2] <<" "<< Evtot[3] \
    <<" "<< Evtot[4] <<" "<< Evtot[5] <<" "<< Evtot[6] << endl;
  //myfile2 << mS<<" "<<mN <<" "<< ctN<<" "<< nEvs2 <<" "<< Evtot2[0] <<" "<< Evtot2[1] <<" "<< Evtot2[2] <<" "<< Evtot2[3] \
    <<" "<< Evtot2[4] <<" "<< Evtot2[5] <<" "<< Evtot2[6] << endl;

  /////########### Get final yield and ratio if excluded
  vector<double> stm{36.7,14.6,15.6,15.1,35.3,20.7,5.61};

  double max = *max_element(Evtot.begin(), Evtot.end());
  int maxid = max_element(Evtot.begin(),Evtot.end()) - Evtot.begin();
  double effi1 = max / nEvs;
  double max2 = *max_element(Evtot2.begin(), Evtot2.end());
  int maxid2 = max_element(Evtot2.begin(),Evtot2.end()) - Evtot2.begin();
  double effi2 = max2 / nEvs2;

  double xs = 1000;
  double lumi = 16.1; //fb-1
  double ny1 = xs * lumi * effi1;
  double ny2 = xs * lumi * effi2;
  
  double exlcms = stm[maxid];
  double exlcms2 = stm[maxid2];
  double rat1 = ny1/exlcms;
  double rat2 = ny2/exlcms2;

  //cout <<"r: "<<endl;
  cout <<"SetNumber:"<< maxid+1 <<" ; signal acceptance: " << effi1<< endl;
  cout <<"CMS, L=16.1fb-1, 95CL excluded: "<< exlcms<< endl;
  cout <<"xsec=1pb, signal events: "<< ny1<<" ; r="<< rat1 << endl;  

  //cout <<"It5: "<<endl; 
  //cout <<"SetNumber:"<< maxid2+1 <<" ; signal acceptance: " << effi2<< endl;
  //cout <<"CMS, L=16.1fb-1, 95CL excluded: "<< exlcms2<< endl;
  //cout <<"xsec=1pb, signal events: "<< ny2<<" ; r="<< rat2 << endl;  

  pythia.stat();

	//myfile1.close();
  //myfile2.close();
  return 0;
}
