#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include<TF1.h>
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
  bool DEBUG=false;
  //counter
  int nHgc=0, nAhc=0, nrechits=0;
  Long64_t jentry;
  //HGCNtupleVariables::init_piTree();
  char* outFileName = new char[1000];
  sprintf(outFileName,"./FullHGCAlSS_ntuple_8k_pion10to350GeV_model2_version73.root");//,energy);
  TFile* outfile = TFile::Open(outFileName,"recreate");
  pion_tree = new TTree("pion_variables_v1","variables for pion analysis v1");
  HGCNtupleVariables::init_piTree();
  int count[350]={};
  float lambda[79];
  nentries=100;
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
      // apparently code is not running beyond this point ///
      pi_Event=0;
      pi_nHadrons[47]={};
      pi_nGammas[47]={};
      pi_nMuons[47]={};
      pi_nElectrons[47]={};
      pi_genId=0.0;
      pi_genEn=0.0;
      pi_genEt=0.0;
      pi_genEta=0.0;
      pi_genPhi=0.0;
      pi_nhits=0;
      pi_hit_sithick[5621]={};
      pi_hit_si[5621]={};
      pi_hit_lay[5621]={};
      pi_hit_x[5621]={};
      pi_hit_y[5621]={};
      pi_hit_z[5621]={};
      pi_hit_men[5621]={};
      pi_hit_en[5621]={};
      pi_hit_endens[5621]={};
      pi_hit_dR[5621]={};
      pi_hit_dRho[5621]={};
      pi_hit_miss[47]={};
      pi_si_clustsumen[5][47]={};
      pi_sci_clustsumen[5][47]={};
      pi_si_sumen[47]={};
      pi_sci_sumen[47]={};
      pi_total_sim_men[47]={};
      pi_avg_emFrac[47]={};
      pi_avg_hadFrac[47]={};

      for (int il=0; il<47;il++)
	{
	  cout<<si_sumen[0]<<"\t"<<si_sumen[1]<<"\t"<<si_sumen[47]<<endl;
	  // pi_si_sumen[il] = si_sumen[il];
	  // pi_sci_sumen[il] = sci_sumen[il];
	}

      h_beamenergy->Fill(genEt);
      h_particle->Fill(genId);
      h_true_beamenergy->Fill(genEn);
      for (int i=10;i<351;i++)
      	{
      	  if(int(genEn)==i)
      	    {
      	      count[i]++;
      	    }
      	}

      // pi_Event=Event;
      
      // pi_nHadrons=nHadrons;
      // pi_nGammas=nGammas;
      // pi_nMuons=nMuons;
      // pi_nElectrons=nElectrons;
      // pi_genId=genId;
      // pi_genEn=genEn;
      // pi_genEt=genEt;
      // pi_genEta=genEta;
      // pi_genPhi=genPhi;
      // pi_nhits=nhits;
      // pi_hit_sithick=hit_sithick;
      // pi_hit_si = hit_si;
      // pi_hit_lay= hit_lay;
      // pi_hit_x = hit_x;
      // pi_hit_y=hit_y;
      // pi_hit_z = hit_z;
      // pi_hit_men = hit_men;
      // pi_hit_en=hit_en;
      // pi_hit_endens = hit_endens;
      // pi_hit_dR = hit_dR;
      // pi_hit_dRho = hit_dRho;
      // pi_hit_miss= hit_miss;
      // pi_si_clustsumen = si_clustsumen;
      // pi_sci_clustsumen = sci_clustsumen;
      // pi_si_sumen= si_sumen;
      // pi_sci_sumen = sci_sumen;
      // pi_total_sim_men = total_sim_men;
      // pi_avg_emFrac = avg_emFrac;
      // pi_avg_hadFrac = avg_hadFrac;

      //      pion_tree->Fill();
      
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
  // char* name1 = new char[1000];
  // sprintf(name1,"./events_vs_energy.txt");
  // std::ofstream file_H_;
  // file_H_.open(name1,ios::out);

  // cout<<"Got Out "<<endl;
  // for (int i=10;i<351;i++)
  //   {
  //     file_H_<<i<<"\t"<<count[i]<<"\n";
  //   }
  // cout<<count_fh<<endl;
  // cout<<count_badtrack<<endl;
  // oFile->cd();
  // oFile->Write();
  // oFile->Close();

}
