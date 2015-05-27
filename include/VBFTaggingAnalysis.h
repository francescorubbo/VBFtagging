//-*-c++-*-

#ifndef  VBFTaggingAnalysis_H
#define  VBFTaggingAnalysis_H

#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "VBFTaggingInfo.h"
#include "Configuration.h"
#include "Definitions.h"

using namespace std;
using namespace fastjet;

typedef vector<float> treeBranch;

class VBFTaggingAnalysis{
 private:
  int  ftest;
  bool  fDebug;
  string fOutName;
  
  TFile *tF;
  TTree *tT;

  //these need to be ptrs because constructor must be called
  unique_ptr<JetDefinition> jetDef;
  unique_ptr<AreaDefinition> active_area;
  unique_ptr<GridMedianBackgroundEstimator> bge;
  unique_ptr<FunctionOfPseudoJet<double> > rescaling;
  unique_ptr<Selector> select_fwd;
    
  // Tree Vars ---------------------------------------
  int fTEventNumber;
  int fTNPV;
  
  treeBranch *jpt;
  treeBranch *jphi;
  treeBranch *jeta;
  treeBranch *jtruth;
  
  treeBranch *truejpt;
  treeBranch *truejphi;
  treeBranch *truejeta;

  Pythia8::Pythia *_pythiaHS;
  Pythia8::Pythia *_pythiaPU;
  
  void DeclareBranches();
  void ResetBranches();
  void FillTree(JetVector jets, JetVector truthJets);
  void FillTruthTree(JetVector jets);
  double TruthFrac(PseudoJet jet, JetVector truthJets);
  bool Ignore(Pythia8::Particle &p);

  //Jet selection functions
  void selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets);
  
 public:
  VBFTaggingAnalysis (Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q);
  ~VBFTaggingAnalysis ();
  
  void AnalyzeEvent(int iEvt, int NPV);
  void Initialize();
  
  //settings (call before initialization)
  void Debug(int debug){fDebug = debug;}
  void SetOutName(string outname){fOutName = outname;}
  
};

#endif

