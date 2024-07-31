#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include<TF1.h>
#define MAXHITS 50000
#define MAXSEC 1000

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include<map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TVector2.h"
#include "TRandom2.h"
#include "TSystem.h"


using namespace std;



// chip 3022,44,3028




int main(int argc, char* argv[])//, int argvv[])
{

  if (argc < 3) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *energy = argv[3];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, energy);
  //cout << "dataset " << data << " " << endl;
  //cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;
  // int chi2_method = atoi(energy);
  hgcOctTB.EventLoop(energy);//, min_, max_);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *energy) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;
  // int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
  // 		    82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
  // 		    154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
  // 		    226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
  // 		    298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  int Elist[8]={20,50,80,100,120,200,250,300};
  cout << "nentries " << nentries << endl;
  //  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;

  int decade = 0;
  char* outFileName = new char[1000];
  int en = atoi(energy);
  sprintf(outFileName,"./SkimmedFiles_wFullHGCAlSS_ntuple_model2_version73_Enegry%d.root",en);                                                
 int pNPhysLayers =47;
 TFile* outfile = TFile::Open(outFileName,"recreate");
 outfile->cd();
  TTree *hitTree=new TTree("hits","hits");
  Int_t pi_genId=0.0;
  Int_t pi_nHadrons[pNPhysLayers]={},pi_nGammas[pNPhysLayers]={},pi_nMuons[pNPhysLayers]={},pi_nElectrons[pNPhysLayers]={};
  Float_t pi_genEn=0.0,pi_genEt=0.0,pi_genEta=0.0,pi_genPhi=0.0;
  Int_t pi_nhits=0, pi_Event=0;
  Bool_t  pi_hit_si[MAXHITS]={};
  Int_t   pi_hit_lay[MAXHITS]={},pi_hit_sithick[MAXHITS]={};
  Float_t pi_hit_x[MAXHITS]={},pi_hit_y[MAXHITS]={},pi_hit_z[MAXHITS]={},pi_hit_men[MAXHITS]={},pi_hit_en[MAXHITS]={},pi_hit_time[MAXHITS]={},pi_hit_endens[MAXHITS]={},pi_hit_dR[MAXHITS]={},pi_hit_dRho[MAXHITS]={};
  Float_t pi_si_sumen[pNPhysLayers]={},pi_sci_sumen[pNPhysLayers]={};
  Float_t pi_si_clustsumen[5][pNPhysLayers];
  Float_t pi_sci_clustsumen[5][pNPhysLayers];
  pi_si_clustsumen[5][pNPhysLayers]={};
  pi_sci_clustsumen[5][pNPhysLayers]={};
  Float_t pi_hit_miss[pNPhysLayers]={};
  Float_t pi_total_sim_men[pNPhysLayers]={},pi_avg_emFrac[pNPhysLayers]={},pi_avg_hadFrac[pNPhysLayers]={};
                                                
  Float_t pi_avg_emFrac_total, pi_avg_hadFrac_total, pi_EneTotal,pi_lay_emEne[pNPhysLayers],pi_lay_hadEne[pNPhysLayers];
  Float_t pi_particle_charge[MAXSEC], pi_particle_kin[MAXSEC], pi_particle_x[MAXSEC], pi_particle_y[MAXSEC], pi_particle_z[MAXSEC];

  // Float_t pi_particle_x,pi_particle_y,pi_particle_z,pi_particle_charge,pi_particle_kin,pi_particle_track_id,pi_particle_parent_id;
  Int_t pi_lay_num;
  Int_t pi_nSec;
  Float_t pi_sec_pdgID[MAXSEC], pi_sec_charge[MAXSEC], pi_sec_kin[MAXSEC], pi_int_x[MAXSEC], pi_int_y[MAXSEC], pi_int_z[MAXSEC];

  Int_t pi_nparticle;
  Float_t pi_particle_pdgID[MAXSEC], pi_article_charge[MAXSEC], pi_article_kin[MAXSEC], pi_article_x[MAXSEC], pi_article_y[MAXSEC], pi_article_z[MAXSEC];
  Float_t pi_particle_process_id[MAXSEC], pi_particle_parent_id[MAXSEC], pi_particle_track_id[MAXSEC];
  std::vector<std::string> *pi_particle_creator_process;
  Float_t pi_measuredE[pNPhysLayers], pi_absorberE[pNPhysLayers], pi_totalE[pNPhysLayers];
  Double_t pi_leadingHadKE;
  Double_t pi_E_genKin;
  Double_t pi_nsc_total_KE;
  Double_t pi_Allsec_pions;
  Double_t pi_Allsec_Kaons;
  Double_t pi_Allseckin;
  int Nrechit_trimAhcal;


  hitTree->Branch("Event", Event,"pi_Event/I");
  hitTree->Branch("nHadrons",    pi_nHadrons,   Form("pi_nHadrons[%d]/I",pNPhysLayers) );
  hitTree->Branch("nGammas",     pi_nGammas,    Form("pi_nGammas[%d]/I",pNPhysLayers) );
  hitTree->Branch("nMuons",      pi_nMuons,     Form("pi_nMuons[%d]/I",pNPhysLayers) );
  hitTree->Branch("nElectrons",  pi_nElectrons, Form("pi_nElectrons[%d]/I",pNPhysLayers) );
  hitTree->Branch("genId",  &pi_genId,   "pi_genId/I");
  hitTree->Branch("genEn",  &pi_genEn,   "pi_genEn/F");
  hitTree->Branch("genEt",  &pi_genEt,   "pi_genEt/F");
  hitTree->Branch("genEta", &pi_genEta,  "pi_genEta/F");
  hitTree->Branch("genPhi", &pi_genPhi,  "pi_genPhi/F");
  hitTree->Branch("nhits",  &pi_nhits,   "pi_nhits/I");
  hitTree->Branch("hit_sithick",  pi_hit_sithick,  "pi_hit_sithick[pi_nhits]/I");
  hitTree->Branch("hit_si",  pi_hit_si,  "pi_hit_si[pi_nhits]/O");
  hitTree->Branch("hit_lay", pi_hit_lay, "pi_hit_lay[pi_nhits]/I");
  hitTree->Branch("hit_x",   pi_hit_x,   "pi_hit_x[pi_nhits]/F");
  hitTree->Branch("hit_y",   pi_hit_y,   "pi_hit_y[pi_nhits]/F");
  hitTree->Branch("hit_z",   pi_hit_z,   "pi_hit_z[pi_nhits]/F");
  hitTree->Branch("hit_men",  pi_hit_men,  "pi_hit_men[pi_nhits]/F");
  hitTree->Branch("hit_en",  pi_hit_en,  "pi_hit_en[pi_nhits]/F");
  hitTree->Branch("hit_time",pi_hit_time,"pi_hit_time[pi_nhits]/F");
  hitTree->Branch("hit_endens",  pi_hit_endens,  "pi_hit_endens[pi_nhits]/F");
  hitTree->Branch("hit_dR",  pi_hit_dR,  "pi_hit_dR[pi_nhits]/F");
  hitTree->Branch("hit_dRho", pi_hit_dRho, "pi_hit_dRho[pi_nhits]/F");
  hitTree->Branch("hit_miss", pi_hit_miss, Form("pi_hit_miss[%d]/F",pNPhysLayers) );
  hitTree->Branch("si_clustsumen", pi_si_clustsumen, Form("pi_si_clustsumen[5][%d]/F",pNPhysLayers) );
  hitTree->Branch("sci_clustsumen", pi_sci_clustsumen, Form("pi_sci_clustsumen[5][%d]/F",pNPhysLayers) );
  hitTree->Branch("si_sumen", pi_si_sumen, Form("pi_si_sumen[%d]/F",pNPhysLayers));
  hitTree->Branch("sci_sumen", pi_sci_sumen, Form("pi_sci_sumen[%d]/F",pNPhysLayers));
  hitTree->Branch("total_sim_men",  pi_total_sim_men,  Form("pi_total_sim_men[%d]/F",pNPhysLayers));
  hitTree->Branch("avg_emFrac",  pi_avg_emFrac,  Form("pi_avg_emFrac[%d]/F",pNPhysLayers));
  hitTree->Branch("avg_hadFrac",  pi_avg_hadFrac,  Form("pi_avg_hadFrac[%d]/F",pNPhysLayers));
  hitTree->Branch("lay_emEne",  pi_lay_emEne,  Form("pi_lay_emEne[%d]/F",pNPhysLayers));
  hitTree->Branch("lay_hadEne",  pi_lay_hadEne,  Form("pi_lay_hadEne[%d]/F",pNPhysLayers));
  hitTree->Branch("avg_emFrac_total",  &pi_avg_emFrac_total,  "pi_avg_emFrac_total/F");
  hitTree->Branch("avg_hadFrac_total",  &pi_avg_hadFrac_total,  "pi_avg_hadFrac_total/F");

  hitTree->Branch("measuredE",  pi_measuredE,  Form("pi_measuredE[%d]/F",pNPhysLayers));
  hitTree->Branch("absorberE",  pi_absorberE,  Form("pi_absorberE[%d]/F",pNPhysLayers));
  hitTree->Branch("totalE",  pi_totalE,  Form("pi_totalE[%d]/F",pNPhysLayers));

  hitTree->Branch("nSec",  &pi_nSec,   "pi_nSec/I");
  hitTree->Branch("int_x",  pi_int_x,  "pi_int_x[pi_nSec]/F");
  hitTree->Branch("int_y",  pi_int_y,   "pi_int_y[pi_nSec]/F");
  hitTree->Branch("int_z",  pi_int_z,   "pi_int_z[pi_nSec]/F");
  hitTree->Branch("sec_pdgID",   pi_sec_pdgID,   "pi_sec_pdgID[pi_nSec]/F");
  hitTree->Branch("sec_charge",   pi_sec_charge,   "pi_sec_charge[pi_nSec]/F");
  hitTree->Branch("sec_kin",   pi_sec_kin,   "pi_sec_kin[pi_nSec]/F");

  hitTree->Branch("nparticle",  &pi_nparticle,   "pi_nparticle/I");
  hitTree->Branch("particle_x",  pi_particle_x,   "pi_particle_x[pi_nparticle]/F");
  hitTree->Branch("particle_y",  pi_particle_y,   "pi_particle_y[pi_nparticle]/F");
  hitTree->Branch("particle_z",  pi_particle_z,   "pi_particle_z[pi_nparticle]/F");
  hitTree->Branch("particle_pdgID",   pi_particle_pdgID,   "pi_particle_pdgID[pi_nparticle]/F");
  hitTree->Branch("particle_charge",   pi_particle_charge,   "pi_particle_charge[pi_nparticle]/F");
  hitTree->Branch("particle_kin",   pi_particle_kin,   "pi_particle_kin[pi_nparticle]/F");
  hitTree->Branch("particle_process_id",   pi_particle_process_id,   "pi_particle_process_id[pi_nparticle]/F");
  hitTree->Branch("particle_parent_id",   pi_particle_parent_id,   "pi_particle_parent_id[pi_nparticle]/F");
  hitTree->Branch("particle_track_id",   pi_particle_track_id,   "pi_particle_track_id[pi_nparticle]/F");
  hitTree->Branch("particle_creator_process",   &pi_particle_creator_process);
  hitTree->Branch("leadingHadKE", &pi_leadingHadKE);
  hitTree->Branch("E_genKin", &pi_E_genKin);
  hitTree->Branch("nsc_total_KE",&pi_nsc_total_KE);
  hitTree->Branch("Allsec_pions",&pi_Allsec_pions);
  hitTree->Branch("Allsec_Kaons",&pi_Allsec_Kaons);
  hitTree->Branch("Allseckin",&pi_Allseckin);
  
  hitTree->SetAutoFlush(1000);
  //  hitTree->SetDirectory(outputFile);

  bool DEBUG=false;
  //counter
  int nHgc=0, nAhc=0, nrechits=0;
  Long64_t jentry;

  int count[350]={};
  bool flag_en[350]={};
  float lambda[79];
  //nentries=100;
  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
  {
    //   for (jentry=0; jentry<10000;jentry++,hgc_jentry++) {
    // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
      if(jentry%10000==0)
	
	cout<<"entre:  "<<jentry<<endl;
      // ===============read this entry == == == == == == == == == == ==

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { break; cout<<"Breaking"<<endl;}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<"insides the event loop: check1"<<endl;

     ////   MESSAGE ////
      int i_en =0;
      //      if(genEn!=200) continue;
      // apparently code is not running beyond this point ///
      // for (int i=10;i<351;i++)
      //   {
      //     if(int(genEn)==i)
      //       {
      //         // if(count[i]>8000)
      //         //   break;

      //         //hitTree->Fill();                                                                                                                 
      // 	      i_en = i;
      //         count[i]++;
      //       }
      //   }
      // if(count[i_en]>8000)
      // 	{//cout<<"final"<<endl;
      // 	  flag_en[i_en] = true;
      // 	  //eak;
      // 	}
      //      cout<<i_en<<"\t"<<count[i_en]<<endl;
      h_beamenergy->Fill(genEn);
      h_particle->Fill(genId);
      //      h_true_beamenergy->Fill(genEn);

      // //initialize
      if(DEBUG)
	cout<<jentry<<" just before the branch initialisation  "<<endl;
       pi_genId=0.0;
       pi_nHadrons[pNPhysLayers]={},pi_nGammas[pNPhysLayers]={},pi_nMuons[pNPhysLayers]={},pi_nElectrons[pNPhysLayers]={};
       pi_genEn=0.0,pi_genEt=0.0,pi_genEta=0.0,pi_genPhi=0.0;
       pi_nhits=0, pi_Event=0;
       pi_hit_si[MAXHITS]={};
       pi_hit_lay[MAXHITS]={},pi_hit_sithick[MAXHITS]={};
       pi_hit_x[MAXHITS]={},pi_hit_y[MAXHITS]={},pi_hit_z[MAXHITS]={},pi_hit_men[MAXHITS]={},pi_hit_en[MAXHITS]={},pi_hit_endens[MAXHITS]={};
       pi_hit_dR[MAXHITS]={},pi_hit_dRho[MAXHITS]={}, pi_hit_time[MAXHITS]={};
       pi_si_sumen[pNPhysLayers]={},pi_sci_sumen[pNPhysLayers]={};
       pi_si_clustsumen[5][pNPhysLayers]={};
       pi_sci_clustsumen[5][pNPhysLayers]={};
       pi_hit_miss[pNPhysLayers]={};
       pi_total_sim_men[pNPhysLayers]={},pi_avg_emFrac[pNPhysLayers]={},pi_avg_hadFrac[pNPhysLayers]={};
       pi_Event=0;
       pi_genEn=0.0,pi_genEt=0.0,pi_genEta=0.0,pi_genPhi=0.0;
       pi_leadingHadKE = 0;
       pi_E_genKin = 0;
       pi_nsc_total_KE = 0;
       pi_Allsec_pions = 0;
       pi_Allsec_Kaons = 0;
       pi_Allseckin =0;


       // additonal variables
       pi_EneTotal=0;
       pi_lay_emEne[pNPhysLayers]={},pi_avg_hadFrac[pNPhysLayers]={};
       pi_avg_emFrac_total=0, pi_avg_hadFrac_total=0,pi_EneTotal=0;
       pi_particle_charge[MAXSEC]={},pi_particle_kin[MAXSEC]={},pi_particle_x[MAXSEC]={},pi_particle_y[MAXSEC]={}, pi_particle_z[MAXSEC]={};
       pi_lay_num=0;
       pi_nSec=0, pi_sec_pdgID[MAXSEC]={},pi_sec_charge[MAXSEC]={},pi_sec_kin[MAXSEC]={};
       pi_int_x[MAXSEC]={},pi_int_y[MAXSEC]={},pi_int_z[MAXSEC]={};
       pi_nparticle=0;
       pi_particle_pdgID[MAXSEC]={},pi_particle_process_id[MAXSEC]={},pi_particle_parent_id[MAXSEC]={},pi_particle_track_id[MAXSEC]={};
       pi_particle_creator_process->clear();
       pi_measuredE[pNPhysLayers]={};
       pi_absorberE[pNPhysLayers]={};
       pi_totalE[pNPhysLayers]={};
      // //filling the branches
       if(DEBUG)
	 cout<<jentry<<" just after finishing intialisation  "<<endl;

      pi_Event= Event;
      pi_genEn = genEn;
      pi_genEt = genEt;
      pi_genEta = genEta;
      pi_genPhi = genPhi;
      pi_genId= genId;
      //additonal variables
      pi_EneTotal = EneTotal;
      pi_avg_emFrac_total =avg_emFrac_total;
      pi_avg_hadFrac_total = avg_hadFrac_total;
      pi_nSec = nSec;
      // pi_int_x = int_x;
      // pi_int_y = int_y; 
      // pi_int_z = int_z; 
      pi_nparticle = nparticle;
      for(int il=0;il<47;il++)
      	{
      	  pi_hit_miss[il] = hit_miss[il];
      	  pi_total_sim_men[il] = total_sim_men[il];
      	  pi_avg_emFrac[il] = avg_emFrac[il];
      	  pi_avg_hadFrac[il] = avg_hadFrac[il];
	  
      	  pi_nHadrons[il] = nHadrons[il];
      	  pi_nGammas[il]= nGammas[il];
      	  pi_nMuons[il]=nMuons[il];
      	  pi_nElectrons[il] = nElectrons[il];
      	  pi_si_sumen[il] = si_sumen[il];
      	  pi_sci_sumen[il]= sci_sumen[il];
	  pi_lay_emEne[il] = lay_emEne[il];
	  pi_avg_hadFrac[il]= avg_hadFrac[il];
	  pi_measuredE[il]= measuredE[il];
	  pi_absorberE[il]=absorberE[il];
	  pi_totalE[il]=totalE[il];
	  
      	}
      for (int is=0; is<nSec;is++)
	{
	  pi_sec_pdgID[is]= sec_pdgID[is];
	  pi_sec_charge[is]=sec_charge[is];
	  pi_sec_kin[is]= sec_kin[is];
	  pi_int_x[is]= int_x[is];
	  pi_int_y[is]= int_y[is];
	  pi_int_z[is]= int_z[is];
	}
      if(DEBUG)
	cout<<jentry<<" just before filling particle process "<<endl;

      for(int ip=0; ip<nparticle;ip++)
      	{
      	  pi_particle_pdgID[ip]=particle_pdgID[ip];
      	  pi_particle_process_id[ip]=particle_process_id[ip];
      	  pi_particle_parent_id[ip]= particle_parent_id[ip];
      	  pi_particle_track_id[ip] = particle_track_id[ip];
	  //      	  pi_particle_creator_process[ip] = particle_creator_process[ip];
      	  pi_particle_x[ip] = particle_x[ip];
      	  pi_particle_y[ip] = particle_y[ip];
      	  pi_particle_z[ip] = particle_z[ip];
      	  pi_particle_kin[ip]= particle_kin[ip];
	  
      	}
      if(DEBUG)
	cout<<jentry<<" just after finishing intialisation  "<<endl;

      for (int j =0; j<5;j++)
      	{
      	  //cout<<pi_si_clustsumen[4][0]<<"\t"<<pi_si_clustsumen[4][0]<<endl;
      	  for(int il=0;il<47;il++)
      	    {
      	      pi_si_clustsumen[j][il] = si_clustsumen[j][il] ;
      	      pi_sci_clustsumen[j][il] =sci_clustsumen[j][il];
      	    }
      	}
      //pi_nhits = nhits;
      int countt =0;
      for (int i =0; i <nhits;i++)
      	{
	  if(hit_dR[i]<0.5 && hit_men[i]>0.5)
	    {
      	  pi_hit_si[i]= hit_si[i];
      	  pi_hit_lay[i]=hit_lay[i];
      	  pi_hit_sithick[i]=hit_sithick[i];
      	  pi_hit_x[i]= hit_x[i];
      	  pi_hit_y[i]= hit_y[i];
      	  pi_hit_z[i]= hit_z[i];
      	  pi_hit_men[i]= hit_men[i];
      	  pi_hit_en[i]= hit_en[i];
      	  pi_hit_endens[i]=hit_endens[i];
      	  pi_hit_dR[i] = hit_dR[i];
      	  pi_hit_dRho[i] = hit_dRho[i] ;
	  pi_hit_time[i]= hit_time[i];
	  countt ++;
	    }
      	}

      if(DEBUG)
	cout<<jentry<<" just before the pi0 values  "<<endl;

      if (countt==0)
      	continue;
      pi_nhits = countt; 
      double nsc_total_KE = 0.0;
      vector<double> KE_list;
      double leadingHadKE = 0.0;
      double E_gen_kin = 0.0;
      KE_list.clear();
      double E_gen_kin_thres = 0.0;
      double Allsec_pions =0.0 ;
      double Allsec_kaons =0.0 ;
      double AllsecKin=0.0;
      double beamEnergy=0.0;
      int counter=0;

      beamEnergy = genEn;
	//hadronic_int++;
	for(int i = 0; i < nSec; i++) {
	  nsc_total_KE += sec_kin[i];
	  if(abs(sec_pdgID[i]) == 211 && sec_charge[i] == -1) KE_list.push_back(sec_kin[i]);
	}
	leadingHadKE = getLeadingKE(KE_list);
	E_gen_kin = nsc_total_KE - leadingHadKE;
	if(leadingHadKE > 0) {
	  if(E_gen_kin > 10 && false) {
	    cout<<jentry<<" NextClosestLayer is \t "<<endl;
	  }
          if(E_gen_kin < E_gen_kin_thres) continue;
	  //          hard_hadronic_int++;

          if(jentry < 20 && false) {
            cout<<jentry<<" : "<<nSec<<" : "<<nsc_total_KE<< " : "<<leadingHadKE<<" : "<<E_gen_kin<<endl;
	  }
	}

	if(E_gen_kin>beamEnergy) continue;
	bool flag=0;
	for(int i =0;i<nparticle;i++)
	  {
	    int j = i-1;
	    if(j<0) j=0;
	    flag=1;
	    if(i!=0)
	      flag = !((particle_x)[i]==(particle_x)[j] and (particle_y)[i]==(particle_y)[j] and (particle_z)[i]==(particle_z)[j] and (particle_pdgID)[i] ==(particle_pdgID)[j]and (particle_parent_id)[i] == (particle_parent_id)[j]);
	    //      flag =1; //neutralizing the previous condition                                                                                                      
	    // if(flag==1){                                                       
	      
	    if((particle_pdgID)[i]==111 and  (particle_process_id)[i]>=121)// and (particle_process_id)[i]<=151)
	      Allsec_pions+=(particle_kin)[i];
	    if((particle_pdgID)[i]==221 and  (particle_process_id)[i]>=121 )//and  (particle_process_id)[i]<=151)
	      Allsec_kaons+=(particle_kin)[i];
	    if(((particle_pdgID)[i]==221 || (particle_pdgID)[i]==111 ) && (particle_process_id)[i]>=121)// &&  (particle_process_id)[i]<=151)//(!((particle_x)[i]==(particle_x)[j] and (particle_y)[i]==(particle_y)[j]) and (particle_z)[i]==(particle_z)[j]))                                                                                          
	  {
	    AllsecKin+=(particle_kin)[i];
	  }
	  
	  }
	if(Allsec_pions>genEn)
	  {counter++;
	    cout<<"Alpana- AllsecKin "<<AllsecKin<<"\t"<<Allsec_kaons<<"\t"<<Allsec_pions<<"\t"<<leadingHadKE<<"\t"<<nsc_total_KE<<" nsc "<<nSec<<"\t"<<beamEnergy<<endl;
	    continue;
	  }
	for(int ien=0;ien<8;ien++)
	  {
	    if(beamEnergy<=Elist[ien]+2 and beamEnergy>=Elist[ien]-2){
	      h_2D_nsc_correlation[ien]->Fill(E_gen_kin,nSec);
	      h_nsc[ien]->Fill(nSec/beamEnergy);
	      h_seckin[ien]->Fill(nsc_total_KE/beamEnergy);
	      h_had_leading_KE[ien]->Fill(leadingHadKE/beamEnergy);
	      h_gen_kin[ien]->Fill(E_gen_kin/beamEnergy);
	      h_2D_nsc_correlation_frac[ien]->Fill(E_gen_kin/beamEnergy,nSec);
	      h_2D_hlead_sec_frac[ien]->Fill(E_gen_kin/beamEnergy,leadingHadKE/beamEnergy);
	      h_Allsec_kin_kaons[ien]->Fill(Allsec_kaons/beamEnergy);
	      h_Allsec_kin_pions[ien]->Fill(Allsec_pions/beamEnergy);
	      h_Allsec_kin[ien]->Fill(AllsecKin/beamEnergy);
	    }
	  }
	pi_leadingHadKE = leadingHadKE;
	pi_E_genKin = E_gen_kin;
	pi_nsc_total_KE = nsc_total_KE;
	pi_Allsec_pions = Allsec_pions;
	pi_Allsec_Kaons = Allsec_kaons;
	pi_Allseckin =AllsecKin;

		hitTree->Fill();      
	h_true_beamenergy->Fill(genEn);  
      //      pi_nHadrons[pNPhysLayers] = nHadrons[pNPhysLayers];
      // for (int i=10;i<351;i++)
      // 	{
      // 	  if(int(genEn)==i)
      // 	    {
      // 	      if(count[i]>8000)
      // 		break;

      // 	      //hitTree->Fill();
      // 	      count[i]++;
      // 	    }
      // 	}
      //cout<<i_en<<"\t"<<count[i_en]<<endl;
      // if(!flag_en[i_en])
      // 	{
      // 	  hitTree->Fill();

      
      // 	  h_true_beamenergy->Fill(genEn);
      // 	}
      //cout<<sizeof(hit_en[jentry])<<endl;
      // int Nrechits[5]={};
      // int Nrechit_en[5][8]={};
      // float total_Energy[5][8]={};
      // double dr_cutValues[5]={0.5,0.1,0.2,0.3,0.5};
      // int total_nrechits[5][8]={};
      // float total_Energy_EE[5][8]={};
      // float total_Energy_FH[5][8]={};
      // float total_Energy_AH[5][8]={};
      // float total_Energy_EE_mips[5][8]={};
      // float total_Energy_FH_mips[5][8]={};
      // float total_Energy_AH_mips[5][8]={};
      // float total_Energy_mips[5][8]={};
      // float rechit_en=0.0;
      // for(int i=0;i<nhits;i++)
      // 	{
      // 	  // if(hit_men[i]>1000)
      // 	  //   cout<<hit_men[i]<<endl;
      // 	  rechit_en+=hit_men[i];
      // 	  //if(hit_lay[i]<30 && hit_lay[i]>31) continue;
      // 	  for(int j=0;j<1;j++)
      // 	    {
      // 	      // if(j==0)
      // 	      // 	{
      // 	      // 	  Nrechits[j]++;
      // 	      // 	  h_rechitEnergy_Mips[j]->Fill(hit_men[i]);
      // 	      // 	  h_rechitEnergy_MeV[j]->Fill(0.001*hit_en[i]);
      // 	      // 	  h_dR[j]->Fill(hit_dR[i]);
      // 	      // 	  h_dR_rechitEnergy[j]->Fill(hit_men[i],hit_dR[i]);
      // 	      // 	  h_dR_rechitEnergy_GeV[j]->Fill(0.001*hit_en[i],hit_dR[i]);
      // 	      // 	  h_rechitX[j]->Fill(hit_x[i]);
      // 	      // 	  h_rechitY[j]->Fill(hit_y[i]);
      // 	      // 	  h_rechitZ[j]->Fill(hit_z[i]);
      // 	      // 	  for (int k=0;k<85;k++)
      // 	      // 	    {
      // 	      // 	      if(Elist[k]-2<=genEt<=Elist[k]+2)
      // 	      // 	  	{
      // 	      // 	  	  total_nrechits[j][k]++;
      // 	      // 	  	  h_bin_rechitEn_Mips[j][k]->Fill(hit_men[i]);
      // 	      // 	  	  h_bin_rechitEn[j][k]->Fill(0.001*hit_en[i]);
      // 	      // 	  	  h_bin_rechitX[j][k]->Fill(hit_x[i]);
      // 	      // 	  	  h_bin_rechitY[j][k]->Fill(hit_y[i]);
      // 	      // 	  	  h_bin_rechitZ[j][k]->Fill(hit_z[i]);
      // 	      // 	  	  total_Energy[j][k]+=(0.001*hit_en[i]);
      // 	      // 	  	}
      // 	      // 	    }
      // 	      // 	}
      // 	      // else
      // 	      // 	{
      // 	      //cout<<hit_x[i]<<endl;
      // 	      hit_x[i] = 0.1*hit_x[i];
      // 	      hit_y[i] = 0.1*hit_y[i];
      // 	      hit_z[i] = 0.1*hit_z[i];
      // 	      //cout<<hit_x[i]<<endl;
      // 	      if(hit_dR[i]>dr_cutValues[j]) continue;
      // 	      else
      // 		{
      // 		  Nrechits[j]++;
      // 		  h_rechitEnergy_Mips[j]->Fill(hit_men[i]);
      // 		  h_rechitEnergy_MeV[j]->Fill(0.001*hit_en[i]);
      // 		  h_dR[j]->Fill(hit_dR[i]);
      // 		  h_dR_rechitEnergy[j]->Fill(hit_men[i],hit_dR[i]);
      // 		  h_dR_rechitEnergy_GeV[j]->Fill(0.001*hit_en[i],hit_dR[i]);
      // 		  h_rechitX[j]->Fill(hit_x[i]);
      // 		  h_rechitY[j]->Fill(hit_y[i]);
      // 		  h_rechitZ[j]->Fill(hit_z[i]);
      // 		  h_rechitXvsY[j]->Fill(hit_x[i],hit_y[i]);
      // 		  h_rechitXvsZ[j]->Fill(hit_z[i],hit_x[i]);
		  

      // 		  for (int k=0;k<8;k++)
      // 		    {
      // 		      if((Elist[k]-2)<=int(genEn) && int(genEn)<=(Elist[k]+2))
      // 			{
      // 			  //  cout<<"genEn"<<genEn<<"\t"<<k<<"\t"<<Elist[k]<<endl;
      // 			  total_nrechits[j][k]++;
      // 			  total_Energy_mips[j][k]+=hit_men[i];
      // 			  h_bin_rechitEn_Mips[j][k]->Fill(hit_men[i]);
      // 			  h_bin_rechitEn[j][k]->Fill(0.001*hit_en[i]);
      // 			  h_bin_rechitX[j][k]->Fill(hit_x[i]);
      // 			  h_bin_rechitY[j][k]->Fill(hit_y[i]);
      // 			  h_bin_rechitZ[j][k]->Fill(hit_z[i]);

      // 			  total_Energy[j][k]+=(0.001*hit_en[i]);
      // 			  if(hit_lay[i]<=26)
      // 			    {
      // 			      h_bin_rechitEn_Mips_EE[j][k]->Fill(hit_men[i]);
      // 			      h_bin_rechitEn_EE[j][k]->Fill(0.001*hit_en[i]);
      // 			      h_bin_rechitX_EE[j][k]->Fill(hit_x[i]);
      // 			      h_bin_rechitY_EE[j][k]->Fill(hit_y[i]);
      // 			      h_bin_rechitZ_EE[j][k]->Fill(hit_z[i]);
      // 			      total_Energy_EE[j][k]+=(0.001*hit_en[i]);
      // 			      total_Energy_EE_mips[j][k]+=hit_men[i];
      // 			      h_bin_rechitXvsY_EE[j][k]->Fill(hit_x[i],hit_y[i]);
      // 			      h_bin_rechitXvsZ_EE[j][k]->Fill(hit_z[i],hit_x[i]);

      // 			    }
      // 			  else if(hit_lay[i]<=38 && hit_lay[i]>26)
      // 			    {
      // 			      h_bin_rechitEn_Mips_FH[j][k]->Fill(hit_men[i]);
      //                         h_bin_rechitEn_FH[j][k]->Fill(0.001*hit_en[i]);
      //                         h_bin_rechitX_FH[j][k]->Fill(hit_x[i]);
      //                         h_bin_rechitY_FH[j][k]->Fill(hit_y[i]);
      //                         h_bin_rechitZ_FH[j][k]->Fill(hit_z[i]);
      // 			      total_Energy_FH[j][k]+=(0.001*hit_en[i]);
      // 			      total_Energy_FH_mips[j][k]+=hit_men[i];
      // 			      h_bin_rechitXvsY_FH[j][k]->Fill(hit_x[i],hit_y[i]);
      //                         h_bin_rechitXvsZ_FH[j][k]->Fill(hit_z[i],hit_x[i]);

      // 			    }
      // 			  else
      // 			    {
      // 			      h_bin_rechitEn_Mips_AH[j][k]->Fill(hit_men[i]);
      //                         h_bin_rechitEn_AH[j][k]->Fill(0.001*hit_en[i]);
      //                         h_bin_rechitX_AH[j][k]->Fill(hit_x[i]);
      //                         h_bin_rechitY_AH[j][k]->Fill(hit_y[i]);
      //                         h_bin_rechitZ_AH[j][k]->Fill(hit_z[i]);
      // 			      total_Energy_AH[j][k]+=(0.001*hit_en[i]);
      // 			      total_Energy_AH_mips[j][k]+=hit_men[i];
      // 			      h_bin_rechitXvsY_AH[j][k]->Fill(hit_x[i],hit_y[i]);
      //                         h_bin_rechitXvsZ_AH[j][k]->Fill(hit_z[i],hit_x[i]);

      // 			    }
      // 			}
      // 		    }
		    
      // 		}
			  
		
      // 	    }
      // 	}
      // //      cout<<rechit_en<<"\t"<<nhits<<endl;
      // //h_nhits->Fill(nhits);
      // for(int j=0;j<1;j++)
      // 	{
      // 	  h_nhits[j]->Fill(Nrechits[j]);
      // 	  for (int k=0;k<8;k++)
      // 	    {
      // 	      if((Elist[k]-2<=genEn&& genEn<=Elist[k]+2))
      // 	  	{
		  
      // 	  	  h_bin_TotalrechitEn[j][k]->Fill( total_Energy[j][k]);
      // 		  h_bin_TotalrechitEn_EE[j][k]->Fill( total_Energy_EE[j][k]);
      // 		  h_bin_TotalrechitEn_FH[j][k]->Fill( total_Energy_FH[j][k]);
      // 		  h_bin_TotalrechitEn_AH[j][k]->Fill( total_Energy_AH[j][k]);
		  
      // 		  float ratio= total_Energy[j][k]/genEn;
      // 		  float ratio_EE= total_Energy_EE[j][k]/genEn;
      // 		  float ratio_FH= total_Energy_FH[j][k]/genEn;
      // 		  float ratio_AH= total_Energy_AH[j][k]/genEn;
      // 		  h_bin_ratiorechitEn[j][k]->Fill(ratio);
      //             h_bin_ratiorechitEn_EE[j][k]->Fill(ratio_EE);
      //             h_bin_ratiorechitEn_FH[j][k]->Fill(ratio_FH);
      //             h_bin_ratiorechitEn_AH[j][k]->Fill(ratio_AH);
      // 		  h_bin_total_Mips[j][k]->Fill( total_Energy_mips[j][k]);
      // 		    h_bin_total_Mips_EE[j][k]->Fill( total_Energy_EE_mips[j][k]);
      // 		    h_bin_total_Mips_FH[j][k]->Fill( total_Energy_FH_mips[j][k]);
      // 		    h_bin_total_Mips_AH[j][k]->Fill( total_Energy_AH_mips[j][k]);
      // 		    h_bin_total_Mips_EE_v1[j][k]->Fill( total_Energy_EE_mips[j][k]);
      // 		    h_bin_total_Mips_FH_v1[j][k]->Fill( total_Energy_FH_mips[j][k]);

      // 	  	  h_bin_nrechits[j][k]->Fill(total_nrechits[j][k]);
		  
      // 	  	}
      // 	    }
      // 	}
      // total_Energy[5][8]={};
      // total_Energy_EE[5][8]={};
      // total_Energy_FH[5][8]={};
      // total_Energy_AH[5][8]={};
      // total_Energy_mips[5][8]={};
      // total_Energy_EE_mips[5][8]={};
      // total_Energy_FH_mips[5][8]={};
      // total_Energy_AH_mips[5][8]={};


      // totalEnergy_inGeV=0;
      // EnergySum_SSinEE=0;
      // EnergySum_SSinFH=0;
      // Esum_rechits_FH=0;
      // Esum_rechits_EE=0;
      // Esum_rechits_AH=0;
      // Esum_rechits_FH_inGeV=0;
      // Esum_rechits_EE_inGeV=0;
      // Esum_rechits_AH_inGeV=0;
      
      // ////////////////////////////////////////////
      // //            HGCAL Part                  //
      // ////////////////////////////////////////////
      // total_rechits=0;
      // rechits_EE=0;
      // rechits_FH=0;
      // rechits_AH=0;
      // if(event ==3482)
      // 	{	//	continue;
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      //if(DEBUG && jentry > 10) break;
      //      gSystem->Exit(0);
      //if(jentry > 2) break;
      //pion_tree->Fill();
    } // loop overentries
  //  char* name = new char[1000];


  // cout<<endl<<endl;


  // ///////////////////////////////////////////////////////////
  // ///////  E N D     O F     E N T R Y     L O O P     //////
  // ///////////////////////////////////////////////////////////
  // //  gSystem->Exit(0);
  char* name1 = new char[1000];
  sprintf(name1,"./events_vs_energy.txt");
  std::ofstream file_H_;
  file_H_.open(name1,ios::out);

  cout<<"Got Out "<<endl;
  for (int i=10;i<351;i++)
    {
      file_H_<<i<<"\t"<<count[i]<<"\n";
    }
  file_H_.close();
  // cout<<count_fh<<endl;
  // cout<<count_badtrack<<endl;
  hitTree->Write();
  outfile->Close();
  // oFile->cd();
  // oFile->Write();
  // oFile->Close();

}
