//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May  7 12:19:41 2022 by ROOT version 6.24/06
// from TTree hits/hits
// found on file: /eos/cms/store/group/dpg_hgcal/comm_hgcal/kalpana/gitV08-08-00-v3/pi-/Skimmed_Files/Geant4_pi-_version73_12.root
//////////////////////////////////////////////////////////

#ifndef HGCNtupleVariables_h
#define HGCNtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
using namespace std;
// Header file for the classes stored in the TTree if any.

class HGCNtupleVariables {
public :
 HGCNtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~HGCNtupleVariables() { }
  /* void    Init(TTree *tree); */
  /* Bool_t  Notify(); */
  Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TTree* pion_tree;
   Int_t           pi_Event;
   Int_t           pi_nHadrons[47];
   Int_t           pi_nGammas[47];
   Int_t           pi_nMuons[47];
   Int_t           pi_nElectrons[47];
   Int_t           pi_genId;
   Float_t         pi_genEn;
   Float_t         pi_genEt;
   Float_t         pi_genEta;
   Float_t         pi_genPhi;
   Int_t           pi_nhits;
   Int_t  pi_hit_sithick[5621];
   Bool_t  pi_hit_si[5621];
   Int_t   pi_hit_lay[5621];
   Float_t         pi_hit_x[5621];
   Float_t         pi_hit_y[5621];
   Float_t         pi_hit_z[5621];
   Float_t         pi_hit_men[5621];
   Float_t         pi_hit_en[5621];
   Float_t         pi_hit_endens[5621];
   Float_t         pi_hit_dR[5621];
   Float_t         pi_hit_dRho[5621];
   Float_t         pi_hit_miss[47];
   Float_t         pi_si_clustsumen[5][47];
   Float_t         pi_sci_clustsumen[5][47];
   Float_t         pi_si_sumen[47];
   Float_t         pi_sci_sumen[47];
   Float_t         pi_total_sim_men[47];
   Float_t         pi_avg_emFrac[47];
   Float_t         pi_avg_hadFrac[47];
   
   //void init_piTree();
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Int_t           nHadrons[47];
   Int_t           nGammas[47];
   Int_t           nMuons[47];
   Int_t           nElectrons[47];
   Int_t           genId;
   Float_t         genEn;
   Float_t         genEt;
   Float_t         genEta;
   Float_t         genPhi;
   Int_t           nhits;
   Int_t           hit_sithick[5621];                                                                                                
   Bool_t          hit_si[5621];                                                                                                     
   Int_t           hit_lay[5621];                                                                                                    
   Float_t         hit_x[5621];                                                                                                      
   Float_t         hit_y[5621];                                                                                                      
   Float_t         hit_z[5621];                                                                                                      
   Float_t         hit_men[5621];                                                                                                    
   Float_t         hit_en[5621];                                                                                                     
   Float_t         hit_endens[5621];                                                                                                 
   Float_t         hit_dR[5621];                                                                                                     
   Float_t         hit_dRho[5621];                                                                                                   
   Float_t         hit_miss[47];
   Float_t         si_clustsumen[5][47];
   Float_t         sci_clustsumen[5][47];
   Float_t         si_sumen[47];
   Float_t         sci_sumen[47];
   Float_t         total_sim_men[47];
   Float_t         avg_emFrac[47];
   Float_t         avg_hadFrac[47];

   /* Int_t           hit_sithick[2689];   //[nhits] */
   /* Bool_t          hit_si[2689];   //[nhits] */
   /* Int_t           hit_lay[2689];   //[nhits] */
   /* Float_t         hit_x[2689];   //[nhits] */
   /* Float_t         hit_y[2689];   //[nhits] */
   /* Float_t         hit_z[2689];   //[nhits] */
   /* Float_t         hit_men[2689];   //[nhits] */
   /* Float_t         hit_en[2689];   //[nhits] */
   /* Float_t         hit_endens[2689];   //[nhits] */
   /* Float_t         hit_dR[2689];   //[nhits] */
   /* Float_t         hit_dRho[2689];   //[nhits] */
   /* Float_t         hit_miss[47]; */
   /* Float_t         si_clustsumen[5][47]; */
   /* Float_t         sci_clustsumen[5][47]; */
   /* Float_t         si_sumen[47]; */
   /* Float_t         sci_sumen[47]; */
   /* Float_t         total_sim_men[47]; */
   /* Float_t         avg_emFrac[47]; */
   /* Float_t         avg_hadFrac[47]; */

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_nHadrons;   //!
   TBranch        *b_nGammas;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_genId;   //!
   TBranch        *b_genEn;   //!
   TBranch        *b_genEt;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_hit_sithick;   //!
   TBranch        *b_hit_si;   //!
   TBranch        *b_hit_lay;   //!
   TBranch        *b_hit_x;   //!
   TBranch        *b_hit_y;   //!
   TBranch        *b_hit_z;   //!
   TBranch        *b_hit_men;   //!
   TBranch        *b_hit_en;   //!
   TBranch        *b_hit_endens;   //!
   TBranch        *b_hit_dR;   //!
   TBranch        *b_hit_dRho;   //!
   TBranch        *b_hit_miss;   //!
   TBranch        *b_si_clustsumen;   //!
   TBranch        *b_sci_clustsumen;   //!
   TBranch        *b_si_sumen;   //!
   TBranch        *b_sci_sumen;   //!
   TBranch        *b_total_sim_men;   //!
   TBranch        *b_avg_emFrac;   //!
   TBranch        *b_avg_hadFrac;   //!

   //   HGCNtupleVariables(TTree);
   //virtual ~HGCNtupleVariables();
   virtual Int_t    Cut(Long64_t entry);
   //   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void init_piTree();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HGCNtupleVariables_cxx
/* HGCNtupleVariables::HGCNtupleVariables(TTree *tree) : fChain(0)  */
/* { */
/* // if parameter tree is not specified (or zero), connect the file */
/* // used to generate this class and read the Tree. */
/*    if (tree == 0) { */
/*       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/dpg_hgcal/comm_hgcal/kalpana/gitV08-08-00-v3/pi-/Skimmed_Files/Geant4_pi-_version73_12.root"); */
/*       if (!f || !f->IsOpen()) { */
/*          f = new TFile("/eos/cms/store/group/dpg_hgcal/comm_hgcal/kalpana/gitV08-08-00-v3/pi-/Skimmed_Files/Geant4_pi-_version73_12.root"); */
/*       } */
/*       f->GetObject("hits",tree); */

/*    } */
/*    Init(tree); */
/* } */

/* HGCNtupleVariables::~HGCNtupleVariables() */
/* { */
/*    if (!fChain) return; */
/*    delete fChain->GetCurrentFile(); */
/* } */

/* Int_t HGCNtupleVariables::GetEntry(Long64_t entry) */
/* { */
/* // Read contents of entry. */
/*    if (!fChain) return 0; */
/*    return fChain->GetEntry(entry); */
/* } */
Long64_t HGCNtupleVariables::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HGCNtupleVariables::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
  nHadrons[47]={};
  nGammas[47]={};
  nMuons[47]={};
  nElectrons[47]={};
  genId=0;
  genEn=0;
  genEt=0;
  genEta=0;
  genPhi=0;
  nhits=0;
  hit_sithick[5621]={};
  hit_si[5621]={};
  hit_lay[5621]={};
  hit_x[5621]={};
  hit_y[5621]={};
  hit_z[5621]={};
  hit_men[5621]={};
  hit_en[5621]={};
  hit_endens[5621]={};
  hit_dR[5621]={};
  hit_dRho[5621]={};
  hit_miss[47]={};
  si_clustsumen[5][47]={};
  sci_clustsumen[5][47]={};
  si_sumen[47]={};
  sci_sumen[47]={};
  total_sim_men[47]={};
  avg_emFrac[47]={};
  avg_hadFrac[47]={};

   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("nHadrons", nHadrons, &b_nHadrons);
   fChain->SetBranchAddress("nGammas", nGammas, &b_nGammas);
   fChain->SetBranchAddress("nMuons", nMuons, &b_nMuons);
   fChain->SetBranchAddress("nElectrons", nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("genId", &genId, &b_genId);
   fChain->SetBranchAddress("genEn", &genEn, &b_genEn);
   fChain->SetBranchAddress("genEt", &genEt, &b_genEt);
   fChain->SetBranchAddress("genEta", &genEta, &b_genEta);
   fChain->SetBranchAddress("genPhi", &genPhi, &b_genPhi);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("hit_sithick", hit_sithick, &b_hit_sithick);
   fChain->SetBranchAddress("hit_si", hit_si, &b_hit_si);
   fChain->SetBranchAddress("hit_lay", hit_lay, &b_hit_lay);
   fChain->SetBranchAddress("hit_x", hit_x, &b_hit_x);
   fChain->SetBranchAddress("hit_y", hit_y, &b_hit_y);
   fChain->SetBranchAddress("hit_z", hit_z, &b_hit_z);
   fChain->SetBranchAddress("hit_men", hit_men, &b_hit_men);
   fChain->SetBranchAddress("hit_en", hit_en, &b_hit_en);
   fChain->SetBranchAddress("hit_endens", hit_endens, &b_hit_endens);
   fChain->SetBranchAddress("hit_dR", hit_dR, &b_hit_dR);
   fChain->SetBranchAddress("hit_dRho", hit_dRho, &b_hit_dRho);
   fChain->SetBranchAddress("hit_miss", hit_miss, &b_hit_miss);
   fChain->SetBranchAddress("si_clustsumen", si_clustsumen, &b_si_clustsumen);
   fChain->SetBranchAddress("sci_clustsumen", sci_clustsumen, &b_sci_clustsumen);
   fChain->SetBranchAddress("si_sumen", si_sumen, &b_si_sumen);
   fChain->SetBranchAddress("sci_sumen", sci_sumen, &b_sci_sumen);
   fChain->SetBranchAddress("total_sim_men", total_sim_men, &b_total_sim_men);
   fChain->SetBranchAddress("avg_emFrac", avg_emFrac, &b_avg_emFrac);
   fChain->SetBranchAddress("avg_hadFrac", avg_hadFrac, &b_avg_hadFrac);
   Notify();
}
void HGCNtupleVariables::init_piTree(){

  pion_tree->Branch("Event", &pi_Event);
  /* pion_tree->Branch("nHadrons", &pi_nHadrons); */
  /* pion_tree->Branch("nGammas", &pi_nGammas); */
  /* pion_tree->Branch("nMuons", &pi_nMuons); */
  /* pion_tree->Branch("nElectrons", &pi_nElectrons); */
  /* pion_tree->Branch("genId",  &pi_genId); */
  /* pion_tree->Branch("genEn", &pi_genEn); */
  /* pion_tree->Branch("genEt",  &pi_genEt); */
  /* pion_tree->Branch("genEta", &pi_genEta); */
  /* pion_tree->Branch("genPhi", &pi_genPhi); */
  /* pion_tree->Branch("nhits", &pi_nhits); */
  /* pion_tree->Branch("hit_sithick", &pi_hit_sithick); */
  /* pion_tree->Branch("hit_si", &pi_hit_si); */
  /* pion_tree->Branch("hit_lay",  &pi_hit_lay); */
  /* pion_tree->Branch("hit_x",  &pi_hit_x); */
  /* pion_tree->Branch("hit_y", &pi_hit_y); */
  /* pion_tree->Branch("hit_z",  &pi_hit_z); */
  /* pion_tree->Branch("hit_men",  &pi_hit_men); */
  /* pion_tree->Branch("hit_en",  &pi_hit_en); */
  /* pion_tree->Branch("hit_endens",  &pi_hit_endens); */
  /* pion_tree->Branch("hit_dR",  &pi_hit_dR); */
  /* pion_tree->Branch("hit_dRho",  &pi_hit_dRho); */
  /* pion_tree->Branch("hit_miss",  &pi_hit_miss); */
  /* pion_tree->Branch("si_clustsumen",  &pi_si_clustsumen); */
  /* pion_tree->Branch("sci_clustsumen",  &pi_sci_clustsumen); */
  /* pion_tree->Branch("si_sumen",  &pi_si_sumen); */
  /* pion_tree->Branch("sci_sumen",  &pi_sci_sumen); */
  /* pion_tree->Branch("total_sim_men",  &pi_total_sim_men); */
  /* pion_tree->Branch("avg_emFrac",  &pi_avg_emFrac); */
  /* pion_tree->Branch("avg_hadFrac", &pi_avg_hadFrac); */

}

Bool_t HGCNtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HGCNtupleVariables::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HGCNtupleVariables::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HGCNtupleVariables_cxx
