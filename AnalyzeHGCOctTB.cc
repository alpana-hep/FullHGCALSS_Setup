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

  sprintf(outFileName,"./SkimmedFiles_8enPoint_wFullHGCAlSS_ntuple_8kevents_model2_version73_13June2023.root");//,energy);                                                
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
      // apparently code is not running beyond this point ///
      for (int i=10;i<351;i++)
        {
          if(int(genEn)==i)
            {
              // if(count[i]>8000)
              //   break;

              //hitTree->Fill();                                                                                                                 
	      i_en = i;
              count[i]++;
            }
        }
      if(count[i_en]>8000)
	{//cout<<"final"<<endl;
	  flag_en[i_en] = true;
	  //eak;
	}
      //      cout<<i_en<<"\t"<<count[i_en]<<endl;
      h_beamenergy->Fill(genEn);
      h_particle->Fill(genId);
      //      h_true_beamenergy->Fill(genEn);

      // //initialize
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
      // //filling the branches
      pi_Event= Event;
      pi_genEn = genEn;
      pi_genEt = genEt;
      pi_genEta = genEta;
      pi_genPhi = genPhi;
      pi_genId= genId;
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
	  
      	}
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

      if (countt==0)
      	continue;
      pi_nhits = countt; 
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
      if(!flag_en[i_en])
	{
	  hitTree->Fill();

      
	  h_true_beamenergy->Fill(genEn);
	}
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
