
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

void copier(TString kind) {

  TFile* file = TFile::Open(kind+".root");
  TTree* originalTree = (TTree*)file->Get("Low_Tree");
  //TFile* ouput = TFile::Open(kind+"_skim_PFlowResidCut.root","RECREATE");
  TFile* ouput = TFile::Open(kind+"_skim.root","RECREATE");
  //TString cuts = "fabs(particle_px[0]) > 0&&fabs(cell_e[0]) > 0&&fabs(track_qoverp[0]) > 0&&PFlowResid<1000&&PFlowResid>-1000";
  TString cuts = "fabs(supercluster_e[0]) > 0";
  TTree* selectedTree = originalTree->CopyTree(cuts);
  selectedTree->Write();

  int Nb = originalTree->GetEntries();
  int Na = selectedTree->GetEntries();

  std::cout << "For " << kind << " file..." << std::endl;
  std::cout << "N events before: " << Nb << std::endl;
  std::cout << "N events after:  " << Na << std::endl;
  std::cout << "efficiency:      " << Na*1./Nb << std::endl;

}

int copyTree() {

  std::vector<TString> files = {"Photons"};

  for(int i = 0; i<files.size(); i++) copier(files[i]);

  exit(0);

}


