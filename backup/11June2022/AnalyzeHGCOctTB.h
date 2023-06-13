#ifndef AnalyzeHGCOctTB_H
#define AnalyzeHGCOctTB_H
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class AnalyzeHGCOctTB : public HGCNtupleVariables{

 public:
  AnalyzeHGCOctTB(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *energy="");
  ~AnalyzeHGCOctTB();
  //Bool_t   FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList);
   Bool_t   FillChain(TChain *chain, const TString &inputFileList); 
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char* );//, const char *, const char *,const char *);
  void     BookHistogram(const char *, const char* );

  std::vector<bool> *noise_flag;
  TFile *oFile;
  const char *conf_;  
  int inEnergy_;
  int event_count[7]={};
  int count=0, count_afterCuts =0;//,count_nflag=0,count_nchannel=0,count_nchanLayer1=0, count_nflag0=0;
  float  totalEnergy_inGeV=0;
  float  EnergySum_SSinEE=0;
  float  EnergySum_SSinFH=0;
  float Esum_rechits_FH=0;
  float Esum_rechits_EE=0;
  float Esum_rechits_AH=0;
  float Esum_rechits_FH_inGeV=0;
  float Esum_rechits_EE_inGeV=0;
  float Esum_rechits_AH_inGeV=0;
  int total_rechits=0, rechits_EE=0, rechits_FH=0, rechits_AH=0;
  
  int count_FH=0,count_fh=0, count_badtrack=0;//,count_adiMask=0;
  int count_FH_afterCuts=0;
  TH1F* h_beamenergy;
  TH1F* h_particle;
  TH1F* h_true_beamenergy;
  TH1F* h_rechitEnergy_Mips[5];
  TH1F* h_rechitEnergy_MeV[5];
  TH1F* h_dR[5];
  TH1F* h_nhits[5];
  TH2F* h_dR_rechitEnergy[5];
  TH2F* h_dR_rechitEnergy_GeV[5];
  TH1F* h_rechitX[5];
  TH1F* h_rechitY[5];
  TH1F* h_rechitZ[5];
  TH2F* h_rechitXvsY[5];
  TH2F* h_rechitXvsZ[5];
  //Energy wise distributions
  TH1F* h_bin_rechitEn[5][85];
  TH1F* h_bin_rechitEn_Mips[5][85];
  
  TH1F* h_bin_rechitX[5][85];
  TH1F* h_bin_rechitY[5][85];
  TH1F* h_bin_rechitZ[5][85];
  TH1F* h_bin_rechitEn_EE[5][85];
  TH1F* h_bin_rechitEn_Mips_EE[5][85];

  TH1F* h_bin_rechitX_EE[5][85];
  TH1F* h_bin_rechitY_EE[5][85];
  TH1F* h_bin_rechitZ_EE[5][85];

  TH1F* h_bin_rechitEn_FH[5][85];
  TH1F* h_bin_rechitEn_Mips_FH[5][85];

  TH1F* h_bin_rechitX_FH[5][85];
  TH1F* h_bin_rechitY_FH[5][85];
  TH1F* h_bin_rechitZ_FH[5][85];

  TH2F* h_bin_rechitXvsY_FH[5][85];
  TH2F* h_bin_rechitXvsZ_FH[5][85];
  
  TH2F* h_bin_rechitXvsY_EE[5][85];
  TH2F* h_bin_rechitXvsZ_EE[5][85];

  TH1F* h_bin_rechitEn_AH[5][85];
  TH1F* h_bin_rechitEn_Mips_AH[5][85];

  TH1F* h_bin_rechitX_AH[5][85];
  TH1F* h_bin_rechitY_AH[5][85];
  TH1F* h_bin_rechitZ_AH[5][85];

  TH2F* h_bin_rechitXvsY_AH[5][85];
  TH2F* h_bin_rechitXvsZ_AH[5][85];

  TH1F* h_bin_nrechits[5][85];
  TH1F* h_bin_TotalrechitEn[5][85];
  TH1F* h_bin_TotalrechitEn_EE[5][85];
  TH1F* h_bin_TotalrechitEn_FH[5][85];
  TH1F* h_bin_TotalrechitEn_AH[5][85];

  TH1F* h_bin_ratiorechitEn[5][85];
  TH1F* h_bin_ratiorechitEn_EE[5][85];
  TH1F* h_bin_ratiorechitEn_FH[5][85];
  TH1F* h_bin_ratiorechitEn_AH[5][85];

  TH1F* h_bin_total_Mips[5][85];
  TH1F* h_bin_total_Mips_EE[5][85];
  TH1F* h_bin_total_Mips_FH[5][85];
  TH1F* h_bin_total_Mips_AH[5][85];
  TH1F* h_bin_total_Mips_EE_v1[5][85];
  TH1F* h_bin_total_Mips_FH_v1[5][85];


};
#endif

#ifdef AnalyzeHGCOctTB_cxx

void AnalyzeHGCOctTB::BookHistogram(const char *outFileName,  const char* energy) {

  //conf_ = conf;
  oFile = new TFile(outFileName, "recreate");
  //cout<<outFileName<<endl;

  oFile->cd();
  int en = atoi(energy);
  float xlow = 0.0;
  float xhigh = en*100*1.5;
  int xbin = 200;//*0.5*en;
  ////////// Book histograms here  //////////////
   const char *dr_cuts[5]= {"Nocut","dR_01","dR_02","dR_03","dR_05"};//,"dPhi_Met","Met_100","Met_250","Met_600","st_300_Met100","pt_st_Met_250","st_300_Met250","nocut"};

   int Elist[85] =  {20, 50 , 80 , 100 , 120 , 200,  250 , 300 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
		     226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  char name[100],title[100];
  char hname[1000],hname1[1000], hname1_2d[1000];
  sprintf(hname,"Et");
  h_beamenergy = new TH1F(hname,hname,400,0,400);
  sprintf(hname,"E");
  h_true_beamenergy=new TH1F(hname,hname,2000,0,2000);
  sprintf(hname,"h_particle");
  h_particle = new TH1F(hname,hname,400,-400,400);
  for(int i=0;i<1;i++)
    {

      sprintf(hname,"h_rechitEnergy_Mips_%s",dr_cuts[i]);
      h_rechitEnergy_Mips[i]= new TH1F(hname,hname,10000,0,1000);
      sprintf(hname,"h_rechitEnergy_MeVs_%s",dr_cuts[i]);
      h_rechitEnergy_MeV[i]= new TH1F(hname,hname,500,0,500);
      sprintf(hname,"h_dR_%s",dr_cuts[i]);
      h_dR[i]= new TH1F(hname,hname,500,0,5.0);
      sprintf(hname,"h_nhits_%s",dr_cuts[i]);
      h_nhits[i]= new TH1F(hname,hname,500,0,6000);
      sprintf(hname,"h_dR_rechitEnergy_%s",dr_cuts[i]);
      h_dR_rechitEnergy[i]= new TH2F(hname,hname,500,0,5000,500,0,5);
      sprintf(hname,"h_dR_rechitEnergy_GeV_%s",dr_cuts[i]);
      h_dR_rechitEnergy_GeV[i]= new TH2F(hname,hname,500,0,500,500,0,5);
      sprintf(hname,"h_rechitX_%s",dr_cuts[i]);
      h_rechitX[i]= new TH1F(hname,hname,600,-300,300);
      sprintf(hname,"h_rechitY_%s",dr_cuts[i]);
      h_rechitY[i]= new TH1F(hname,hname,600,-300,300);
      sprintf(hname,"h_rechitZ_%s",dr_cuts[i]);
      h_rechitZ[i]= new TH1F(hname,hname,400,200,600);
      sprintf(hname,"h_rechitXvsY_%s",dr_cuts[i]);
      h_rechitXvsY[i]= new TH2F(hname,hname,600,-300,300,600,-300,300);
      sprintf(hname,"h_rechitXvsZ_%s",dr_cuts[i]);
      h_rechitXvsZ[i]= new TH2F(hname,hname,400,200,600,600,-300,300);


      for(int j=0;j<8;j++)
	{
	  sprintf(hname1,"h_bin_rechitEnergy_Mips_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitEn_Mips[i][j]= new TH1F(hname1,hname1,10000,0,1000);
	  sprintf(hname1,"h_bin_rechitEnergy_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitEn[i][j]=new TH1F(hname1,hname1,500,0,500);
	  sprintf(hname1,"h_bin_rechitX_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitX[i][j]= new TH1F(hname1,hname1,600,-300,300);
	  sprintf(hname1,"h_bin_rechitY_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitY[i][j]= new TH1F(hname1,hname1,600,-300,300);
	  sprintf(hname1,"h_bin_rechitZ_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitZ[i][j]= new TH1F(hname1,hname1,400,200,600);

	  sprintf(hname1,"h_bin_rechitX_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitX_EE[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitY_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitY_EE[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitZ_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitZ_EE[i][j]= new TH1F(hname1,hname1,400,200,600);
	  
	  sprintf(hname1,"h_bin_rechitEnergy_Mips_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_Mips_EE[i][j]= new TH1F(hname1,hname1,10000,0,1000);
          sprintf(hname1,"h_bin_rechitEnergy_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_EE[i][j]=new TH1F(hname1,hname1,500,0,500);


          sprintf(hname1,"h_bin_rechitX_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitX_FH[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitY_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitY_FH[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitZ_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitZ_FH[i][j]= new TH1F(hname1,hname1,400,200,600);

          sprintf(hname1,"h_bin_rechitEnergy_Mips_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_Mips_FH[i][j]= new TH1F(hname1,hname1,10000,0,1000);
          sprintf(hname1,"h_bin_rechitEnergy_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_FH[i][j]=new TH1F(hname1,hname1,500,0,500);
	  sprintf(hname,"h_bin_rechitXvsY_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitXvsY_FH[i][j]= new TH2F(hname,hname,600,-300,300,600,-300,300);
          sprintf(hname,"h_bin_rechitXvsZ_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitXvsZ_FH[i][j]= new TH2F(hname,hname,400,200,600,600,-300,300);

	  sprintf(hname1,"h_bin_rechitX_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitX_AH[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitY_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitY_AH[i][j]= new TH1F(hname1,hname1,600,-300,300);
          sprintf(hname1,"h_bin_rechitZ_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitZ_AH[i][j]= new TH1F(hname1,hname1,400,200,600);

          sprintf(hname1,"h_bin_rechitEnergy_Mips_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_Mips_AH[i][j]= new TH1F(hname1,hname1,10000,0,1000);
          sprintf(hname1,"h_bin_rechitEnergy_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitEn_AH[i][j]=new TH1F(hname1,hname1,500,0,500);
          sprintf(hname,"h_bin_rechitXvsY_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitXvsY_AH[i][j]= new TH2F(hname,hname,600,-300,300,600,-300,300);
          sprintf(hname,"h_bin_rechitXvsZ_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_rechitXvsZ_AH[i][j]= new TH2F(hname,hname,400,200,600,600,-300,300);


	  sprintf(hname,"h_bin_rechitXvsY_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitXvsY_EE[i][j]= new TH2F(hname,hname,600,-300,300,600,-300,300);
	  sprintf(hname,"h_bin_rechitXvsZ_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_rechitXvsZ_EE[i][j]= new TH2F(hname,hname,400,200,600,600,-300,300);

	  
	  sprintf(hname1,"h_bin_nrechits_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_nrechits[i][j]= new TH1F(hname1,hname1,600,0,6000);
	  sprintf(hname1,"h_bin_TotalrechitEn_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  
	  int ibin=0,xhigh=0;
	  if(Elist[j]<100)
	    {
	      ibin= 5*Elist[j];
	      xhigh = ibin;
	      ibin=500;
	    }
	  else
	    {
	      ibin = 5.0*Elist[j];
	      xhigh = ibin;
	      ibin=500;
	    }
	  
	  h_bin_TotalrechitEn[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
	  sprintf(hname1,"h_bin_TotalrechitEn_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_TotalrechitEn_EE[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
	  sprintf(hname1,"h_bin_TotalrechitEn_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_TotalrechitEn_FH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
	  sprintf(hname1,"h_bin_TotalrechitEn_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_TotalrechitEn_AH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);

	  ibin = 1000*Elist[j];
	  xhigh = ibin;
	  ibin= 10000;
	  xhigh = 10000;
	  sprintf(hname1,"h_bin_total_Mips_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_total_Mips[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_total_Mips_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_total_Mips_EE[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_total_Mips_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_total_Mips_FH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_total_Mips_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_total_Mips_AH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);

	  /* sprintf(hname1,"h_bin_total_Mips_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]); */
          /* h_bin_total_Mips_EE[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh); */
	  ibin= 40000;
          xhigh = 40000;

	  sprintf(hname1,"h_bin_total_Mips_EE_v1_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_total_Mips_EE_v1[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
	  sprintf(hname1,"h_bin_total_Mips_FH_v1_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_total_Mips_FH_v1[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);

	  ibin = 500;
	  xhigh = 50;
	  sprintf(hname1,"h_bin_ratiorechitEn_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
	  h_bin_ratiorechitEn[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_ratiorechitEn_EE_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_ratiorechitEn_EE[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_ratiorechitEn_FH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_ratiorechitEn_FH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
          sprintf(hname1,"h_bin_ratiorechitEn_AH_%s_TrueEn_%d",dr_cuts[i],Elist[j]);
          h_bin_ratiorechitEn_AH[i][j]= new TH1F(hname1,hname1,ibin,0,xhigh);
	  
	  


	}

    }

}
AnalyzeHGCOctTB::AnalyzeHGCOctTB(const TString &inputFileList, const char *outFileName, const char* energy) {
  //,const char* min_, const char* max_) {

  TChain *tree = new TChain("hits");
  // TChain *tree2 = new TChain("trackimpactntupler/impactPoints");
  // TChain *tree3 = new TChain("bigtree");


  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << std::endl;
  }

  HGCNtupleVariables::Init(tree);//, tree2, tree3);
  //  HGCNtupleVariables::init_piTree();
  BookHistogram(outFileName,energy);
  //  EventLoop(energy);
  // moduleMap_init(config);
  // Alignment_Map_Init();
  // Noise_Map_Init();
  // layerPosition_init();
  // Chi2_Weight_Map_Init(chi2_method);  
}


Bool_t AnalyzeHGCOctTB::FillChain(TChain *chain,  const TString &inputFileList) {//TChain *chain2, TChain *chain3, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
    //if(!infile.good()) break;
    // chain2->Add(buffer.c_str());
    // chain3->Add(buffer.c_str());

  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain3 : " << chain3->GetEntries() << std::endl;

  return kTRUE;
}

Long64_t AnalyzeHGCOctTB::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
      Notify();
  }

  // if (!fChain2) return -5;
  // Long64_t centry2 = fChain2->LoadTree(entry);
  // if (centry2 < 0) return centry2;
  // if (!fChain2->InheritsFrom(TChain::Class()))  return centry2;
  // TChain *chain2 = (TChain*)fChain2;
  // if (chain2->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }

  // if (!fChain3) return -5;
  // Long64_t centry3 = fChain3->LoadTree(entry);
  // if (centry3 < 0) return centry3;
  // if (!fChain3->InheritsFrom(TChain::Class()))  return centry3;
  // TChain *chain3 = (TChain*)fChain3;
  // if (chain3->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }
  
  //if (centry==centry2)
  return centry;
  // cout<<"centry = "<<centry<<endl;
  // if(centry>0)
  //   return centry;
  // else return -1;
}

AnalyzeHGCOctTB::~AnalyzeHGCOctTB() { 

  // if (!fChain || !fChain2) return;
  // delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();

  /* cout<<"in the before fchain"<<endl; */
  if (!fChain) return;
  delete fChain->GetCurrentFile();

  oFile->cd();
  oFile->Write();		/*  */
  oFile->Close();
  /* cout<<"in the before fchain"<<endl; */
}

#endif

/*  LocalWords:  Nrechit EE R1 FH GetXaxis SetTitle Sumw2 TH2F reg3 NRechits
 */
/*  LocalWords:  GetYaxis SetTitleOffset
 */
