//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  6 08:41:43 2025 by ROOT version 6.34.02
// from TTree SubRun/SubRun TTree
// found on file: new_numi_flux_run1_fhc_pandora_ntuple.root
//////////////////////////////////////////////////////////

#ifndef SubRun_h
#define SubRun_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class SubRun {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subRun;
   Float_t         pot;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subRun;   //!
   TBranch        *b_pot;   //!

   SubRun(TTree *tree=0);
   virtual ~SubRun();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SubRun_cxx
SubRun::SubRun(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//
 
   TString MyFile= "RunAll_fhc_new_numi_flux_fhc_pandora_ntuple.root"; 
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(MyFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(MyFile);
      }
      TDirectory * dir = (TDirectory*)f->Get(MyFile+":/nuselection");
      dir->GetObject("SubRun",tree);

   }
   Init(tree);
}

SubRun::~SubRun()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SubRun::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SubRun::LoadTree(Long64_t entry)
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

void SubRun::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subRun", &subRun, &b_subRun);
   fChain->SetBranchAddress("pot", &pot, &b_pot);
   Notify();
}

bool SubRun::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void SubRun::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SubRun::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SubRun_cxx
