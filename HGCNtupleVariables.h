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
  double getLeadingKE(vector<double> list);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Int_t           hit_sithick[5621];   //[nhits]                                                                                                 
   Bool_t          hit_si[5621];   //[nhits]                                                                                                      
   Int_t           hit_lay[5621];   //[nhits]                                                                                                     
   Float_t         hit_x[5621];   //[nhits]                                                                                                       
   Float_t         hit_y[5621];   //[nhits]                                                                                                      
   Float_t EneTotal;
   Float_t avg_emFrac_total;
   Float_t avg_hadFrac_total;

 
   Float_t         hit_z[5621];   //[nhits]                                                                                                       
   Float_t         hit_men[5621];   //[nhits]                                                                                                     
   Float_t         hit_en[5621];   //[nhits]                                                                                                      
   Float_t         hit_endens[5621];   //[nhits]                                                                                                  
   Float_t         hit_dR[5621];   //[nhits]                                                                                                     
Float_t         lay_emEne[47];
   Float_t         lay_hadEne[47];
   //   Float_t         totalE[47];
   Float_t         hit_dRho[5621];   //[nhits]                                                                                                   
 Float_t         hit_time[5621];
   Float_t         hit_miss[47];
   Float_t         si_clustsumen[5][47];
   Float_t         sci_clustsumen[5][47];
   Float_t         si_sumen[47];
   Float_t         sci_sumen[47];
   Float_t         total_sim_men[47];
   Float_t         avg_emFrac[47];
   Float_t         avg_hadFrac[47];
   Float_t         measuredE[47];
   Float_t         absorberE[47];
   Float_t         totalE[47];
   Int_t           nSec;
   Float_t         int_x[999];
   Float_t         int_y[999];
   Float_t         int_z[999];
   Float_t         sec_pdgID[999];   //[nSec]                                                                                                    
   Float_t lay_num;
   Float_t         sec_charge[999];   //[nSec]                                                                                                    
   Float_t         sec_kin[999];   //[nSec]                  
   Int_t           nparticle;
   Float_t         particle_x[999];   //[nparticle]                                                                                               
   Float_t         particle_y[999];   //[nparticle]                                                                                               
   Float_t         particle_z[999];   //[nparticle]                                                                                               
   Float_t         particle_pdgID[999];   //[nparticle]                                                                                           
   Float_t         particle_charge[999];   //[nparticle]                                                                                          
   Float_t         particle_kin[999];   //[nparticle]                                                                                             
   Float_t         particle_process_id[999];   //[nparticle]                                                                                      
   Float_t         particle_parent_id[999];   //[nparticle]                         
   Float_t         particle_track_id[999];   //[nparticle]                                                                                        
   vector<string>  *particle_creator_process;


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
   TBranch        *b_hit_time;
   TBranch        *b_avg_hadFrac;   //!
   TBranch        *b_lay_emEne;   //!                                                                                                             
   TBranch        *b_lay_hadEne;   //!                                                                                                            
   TBranch        *b_avg_emFrac_total;   //!                                                                                                      
   TBranch        *b_avg_hadFrac_total;   //!                                                                                                     
   TBranch        *b_measuredE;   //!                            
   TBranch        *b_absorberE;   //!                                                                                                             
   TBranch        *b_totalE;   //!                                                                                                                
   TBranch        *b_nSec;   //!                                                                                                                  
   TBranch        *b_int_x;   //!                                                                                                                 
   TBranch        *b_int_y;   //!                                                                                                                 
   TBranch        *b_int_z;   //!                                                                                                                 
   TBranch        *b_sec_pdgID;   //!                                                                                                             
   TBranch        *b_sec_charge;   //!                                                                                                            
   TBranch        *b_sec_kin;   //!                
   TBranch        *b_nparticle;   //!                                                                                                             
   TBranch        *b_particle_x;   //!                                                                                                            
   TBranch        *b_particle_y;   //!                                                                                                            
   TBranch        *b_particle_z;   //!                                                                                                            
   TBranch        *b_particle_pdgID;   //!                                                                                                        
   TBranch        *b_particle_charge;   //!                                                                                                       
   TBranch        *b_particle_kin;   //!                                                                                                          
   TBranch        *b_particle_process_id;   //!                                                                                                   
   TBranch        *b_particle_parent_id;   //!                                                                                                    
   TBranch        *b_particle_track_id;   //!                                                                                                     
   TBranch        *b_particle_creator_process;   //!                                                                                              


   //   HGCNtupleVariables(TTree);
   //virtual ~HGCNtupleVariables();
   virtual Int_t    Cut(Long64_t entry);
   //   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   /* virtual void     Loop(); */
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
  /* nHadrons[47]={}; */
  /* nGammas[47]={}; */
  /* nMuons[47]={}; */
  /* nElectrons[47]={}; */
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
  hit_time[5621]={};

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
  EneTotal=0;
  lay_emEne[47]={},avg_hadFrac[47]={};
  avg_emFrac_total=0, avg_hadFrac_total=0,EneTotal=0;
  particle_charge[999]={},particle_kin[999]={},particle_x[999]={},particle_y[999]={}, particle_z[999]={};
  lay_num=0;
  nSec=0, sec_pdgID[999]={},sec_charge[999]={},sec_kin[999]={};
  int_x[999]={},int_y[999]={},int_z[99]={};
  nparticle=0;
  particle_pdgID[999]={},particle_process_id[999]={},particle_parent_id[999]={},particle_track_id[999]={};
  particle_creator_process=0;
  measuredE[47]={};
  absorberE[47]={};
  totalE[47]={};

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
   fChain->SetBranchAddress("hit_time", hit_time, &b_hit_time);
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
   fChain->SetBranchAddress("lay_emEne", lay_emEne, &b_lay_emEne);
   fChain->SetBranchAddress("lay_hadEne", lay_hadEne, &b_lay_hadEne);
   fChain->SetBranchAddress("avg_emFrac_total", &avg_emFrac_total, &b_avg_emFrac_total);
   fChain->SetBranchAddress("avg_hadFrac_total", &avg_hadFrac_total, &b_avg_hadFrac_total);
   fChain->SetBranchAddress("measuredE", measuredE, &b_measuredE);
   fChain->SetBranchAddress("absorberE", absorberE, &b_absorberE);
   fChain->SetBranchAddress("totalE", totalE, &b_totalE);
   fChain->SetBranchAddress("nSec", &nSec, &b_nSec);
   fChain->SetBranchAddress("int_x", &int_x, &b_int_x);
   fChain->SetBranchAddress("int_y", &int_y, &b_int_y);
   fChain->SetBranchAddress("int_z", &int_z, &b_int_z);
   fChain->SetBranchAddress("sec_pdgID", sec_pdgID, &b_sec_pdgID);
   fChain->SetBranchAddress("sec_charge", sec_charge, &b_sec_charge);
   fChain->SetBranchAddress("sec_kin", sec_kin, &b_sec_kin);
   fChain->SetBranchAddress("nparticle", &nparticle, &b_nparticle);
   fChain->SetBranchAddress("particle_x", particle_x, &b_particle_x);
   fChain->SetBranchAddress("particle_y", particle_y, &b_particle_y);
   fChain->SetBranchAddress("particle_z", particle_z, &b_particle_z);
   fChain->SetBranchAddress("particle_pdgID", particle_pdgID, &b_particle_pdgID);
   fChain->SetBranchAddress("particle_charge", particle_charge, &b_particle_charge);
   fChain->SetBranchAddress("particle_kin", particle_kin, &b_particle_kin);
   fChain->SetBranchAddress("particle_process_id", particle_process_id, &b_particle_process_id);
   fChain->SetBranchAddress("particle_parent_id", particle_parent_id, &b_particle_parent_id);
   fChain->SetBranchAddress("particle_track_id", particle_track_id, &b_particle_track_id);
   fChain->SetBranchAddress("particle_creator_process", &particle_creator_process, &b_particle_creator_process);
   Notify();
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
