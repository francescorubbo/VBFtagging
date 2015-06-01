#include "VBFTaggingAnalysis.h"
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

using namespace std;

// Constructor 
VBFTaggingAnalysis::VBFTaggingAnalysis(Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q){

  if((pythiaHS != NULL) and (pythiaPU != NULL)){
    _pythiaHS=pythiaHS;
    _pythiaPU=pythiaPU;
  }
  else{
    cerr << "Invalid Pythia pointer passed to VBFTaggingAnalysis" << endl;
    exit(1);
  }

  fDebug=q.fDebug;
  if(fDebug) 
    cout << "VBFTaggingAnalysis::VBFTaggingAnalysis Start " << endl;

  ftest = 0;
  fOutName = q.outName;

  if(fDebug) 
    cout << "VBFTaggingAnalysis::VBFTaggingAnalysis End " << endl;
}

// Destructor 
VBFTaggingAnalysis::~VBFTaggingAnalysis(){
  if(tT != NULL){
    tT->Write();
    tF->Close();
    delete tF;
  }

  if(jpt != NULL){
    delete jpt;
    delete jphi;
    delete jeta;
    delete jtruth;
    delete fjvt;
    delete truejpt;
    delete truejphi;
    delete truejeta;
  }
}

// Begin method
void VBFTaggingAnalysis::Initialize(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("tree", "Event Tree for VBFTagging");

   const double maxEta = 4.5;
   const double R=0.4;
   const double grid_spacing(0.6);
   
   //suppress fastjet banner
   fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);
   
   jetDef.reset(new JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best));
   active_area.reset(new AreaDefinition(fastjet::active_area));
   bge.reset(new GridMedianBackgroundEstimator(maxEta, grid_spacing));
   rescaling.reset(new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05)); //function for rapidity rescaling of rho
   bge->set_rescaling_class(rescaling.get());
   select_fwd.reset(new Selector(SelectorAbsRapRange(0.,maxEta)));
 
   DeclareBranches();
   
   jpt =      new treeBranch();  
   jphi =     new treeBranch();  
   jeta =     new treeBranch();  
   jtruth =   new treeBranch();
   fjvt =     new treeBranch();

   truejpt =  new treeBranch();  
   truejphi = new treeBranch();  
   truejeta = new treeBranch();  

   ResetBranches();
   
   return;
}

// Analyze
void VBFTaggingAnalysis::AnalyzeEvent(int ievt, int NPV){

  if(fDebug) 
    cout << "VBFTaggingAnalysis::AnalyzeEvent Begin " << endl;
  
  // -------------------------
  if (!_pythiaHS->next()) return;
  if(fDebug) 
    cout << "VBFTaggingAnalysis::AnalyzeEvent Event Number " << ievt << endl;
  
  // reset branches 
  ResetBranches();
  
  // new event-----------------------
  fTEventNumber = ievt;
  JetVector particlesForJets;
  JetVector particlesForJets_np;
  
  //Pileup Loop
  
  fTNPV = NPV;
  
  //Loop over Pileup Events
  for (int iPU = 0; iPU <= NPV; ++iPU) {
        
    //Loop over pileup particles
    for (int i = 0; i < _pythiaPU->event.size(); ++i) {

      if(Ignore(_pythiaPU->event[i]))
	continue;
      
      //Instantiate new pseudojet
      PseudoJet p(_pythiaPU->event[i].px(), 
		  _pythiaPU->event[i].py(), 
		  _pythiaPU->event[i].pz(),
		  _pythiaPU->event[i].e() ); 
                  
      p.set_user_info(new VBFTaggingInfo(_pythiaPU->event[i].id(),
					 i,iPU,true)); 
      particlesForJets.push_back(p); 
    }
    if (!_pythiaPU->next()) continue;
  }
  
  // Particle loop -----------------------------------------------------------
  for (int ip=0; ip<_pythiaHS->event.size(); ++ip){
    
    if(Ignore(_pythiaHS->event[ip]))
      continue;
    
    fastjet::PseudoJet p(_pythiaHS->event[ip].px(), 
			 _pythiaHS->event[ip].py(), 
			 _pythiaHS->event[ip].pz(),
			 _pythiaHS->event[ip].e() ); 
    
    p.set_user_info(new VBFTaggingInfo(_pythiaHS->event[ip].id(),
				       ip,0, false));  
    
    particlesForJets.push_back(p);
    particlesForJets_np.push_back(p);
    
  } // end particle loop -----------------------------------------------

  JetVector selectedJets,selectedTruthJets;
  fastjet::ClusterSequenceArea clustSeq(particlesForJets, *jetDef, *active_area);
  selectJets(particlesForJets,clustSeq,selectedJets);

  fastjet::ClusterSequenceArea clustSeqTruth(particlesForJets_np, *jetDef, *active_area);
  selectedTruthJets = sorted_by_pt(clustSeqTruth.inclusive_jets(10.));
  
  FillTree(selectedJets,selectedTruthJets);
  FillTruthTree(selectedTruthJets);
  
  tT->Fill();
  
  if(fDebug) 
    cout << "VBFTaggingAnalysis::AnalyzeEvent End " << endl;
  
  return;
}

void VBFTaggingAnalysis::selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets){
  try{
    bge->set_particles(particlesForJets);

    fastjet::Subtractor subtractor(bge.get());    
    JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
    JetVector subtractedJets = subtractor(inclusiveJets);

    JetVector allSelectedJets;
    allSelectedJets.clear();
    // allSelectedJets = (*select_fwd)(subtractedJets);
    allSelectedJets = (*select_fwd)(inclusiveJets);
    
    //select jets with pt > 20
    selectedJets.clear();
    for( auto ijet = allSelectedJets.begin(); ijet != allSelectedJets.end(); ++ijet){
      if(ijet->pt() >= 10)
	selectedJets.push_back(*ijet);
    }
  }
  catch(...){
    cerr << "Fastjet error caught in selectJets" << endl;
    exit(20);
  }
}

// worker function to actually perform an analysis
void VBFTaggingAnalysis::FillTree(JetVector jets, JetVector TruthJets){

  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    jpt->push_back(jets[ijet].pt());
    jeta->push_back(jets[ijet].eta());
    jphi->push_back(jets[ijet].phi());
    jtruth->push_back(TruthFrac(jets[ijet],TruthJets));
  }
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    // cout << "fwd jet pt: " << jets[ijet].pt() << " eta: " << jets[ijet].eta() << " phi: "<< jets[ijet].phi() 
    // 	 << " truth: " << jtruth->at(ijet) << endl;
    fjvt->push_back(ComputeFJVT(jets[ijet],jets));
  }

}

void VBFTaggingAnalysis::FillTruthTree(JetVector jets){  
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    truejpt->push_back(jets[ijet].pt());
    truejeta->push_back(jets[ijet].eta());
    truejphi->push_back(jets[ijet].phi());
  }
}

double VBFTaggingAnalysis::TruthFrac(PseudoJet jet, JetVector truthJets){

  double maxFrac=0;

  //for each truth jet
  for (unsigned int tj = 0; tj < truthJets.size(); tj++){
    double ptTruthTot=0;
    // double ptTot=truthJets[tj].pt();
    double ptTot = 0;
    //for each truth particle
    for (unsigned int ti=0; ti < truthJets[tj].constituents().size(); ti++){
      auto truthInfo = truthJets[tj].constituents()[ti].user_info<VBFTaggingInfo>();
      int truthID = truthInfo.pythia_id();
      ptTot += truthJets[tj].constituents()[ti].pt();
      
      //for each particle in the main jet
      for (unsigned int i=0; i < jet.constituents().size(); i++){
	auto pInfo = jet.constituents()[i].user_info<VBFTaggingInfo>();
	// if it is the identical particle, only use if we have timing info
	if ((not pInfo.pileup()) and (truthID == pInfo.pythia_id())){
	  ptTruthTot += truthJets[tj].constituents()[ti].pt();
	}
      }
    }
    
    double frac=ptTruthTot/ptTot;
    if(frac > maxFrac)
      maxFrac=frac;
  }
  
  return maxFrac;
}

bool VBFTaggingAnalysis::Ignore(Pythia8::Particle &p){
  if (!p.isFinal() )      
    return true;
  switch(abs(p.id())){
  case 12:
  case 13:
  case 14:
  case 16:
    return true;
  default:
    return false;
  }
}

// declare branches
void VBFTaggingAnalysis::DeclareBranches(){
   // Event Properties 
  gROOT->ProcessLine("#include <vector>");
  
  tT->Branch("EventNumber",&fTEventNumber,"EventNumber/I");
  tT->Branch("NPV",&fTNPV,"NPV/I");
  tT->Branch("jpt","std::vector<float>",&jpt);
  tT->Branch("jphi","std::vector<float>",&jphi);
  tT->Branch("jeta","std::vector<float>",&jeta);
  tT->Branch("jtruth","std::vector<float>",&jtruth);
  tT->Branch("fjvt","std::vector<float>",&fjvt);

  tT->Branch("truejpt", "std::vector<float>",&truejpt);
  tT->Branch("truejphi","std::vector<float>",&truejphi);
  tT->Branch("truejeta","std::vector<float>",&truejeta);
  
  return;
}

// resets vars
void VBFTaggingAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV = -1;

      jpt->clear();
      jphi->clear();
      jeta->clear();
      jtruth->clear();
      fjvt->clear();
      truejpt->clear();
      truejphi->clear();
      truejeta->clear();
}

double VBFTaggingAnalysis::ComputeFJVT(PseudoJet jet, JetVector jets){

  double fjvt = -1.;
  if(abs(jet.eta())<2.5) return fjvt;

  normal gaus(jet.pt(),0.2*jet.pt());
  double totgaus = 0.;
  //for each truth jet
  for (unsigned int ijet = 0; ijet < jets.size(); ijet++){
    if(jtruth->at(ijet)>0.5) continue;
    if(abs(jets[ijet].eta())>2.5) continue;
    PseudoJet b2bjet;
    b2bjet.reset_PtYPhiM(jet.pt(),jet.eta(),-1*jet.phi(),jet.m());
    double phidist = abs(jets[ijet].delta_phi_to(b2bjet));
    if(phidist>0.2) continue;
    totgaus += pdf(gaus,jets[ijet].pt());
    // cout << "ctl jet pt: " << jets[ijet].pt() << " eta: " << jets[ijet].eta() << " phi: " << phidist
    // 	 << " gaus " << pdf(gaus,jets[ijet].pt()) << endl;
  }
  if(totgaus>0.) fjvt = totgaus;

  return fjvt;
}
