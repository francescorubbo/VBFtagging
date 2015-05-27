#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "Definitions.h"
#include "VBFTaggingAnalysis.h"
#include "Configuration.h"

using namespace std;

int main(int argc, char* argv[]){

  //read options, configure Pythia
  Configuration settings(argc, argv);
  std::shared_ptr<Pythia8::Pythia> pythiaHS(new Pythia8::Pythia("../xmldoc",false));
  std::shared_ptr<Pythia8::Pythia> pythiaPU(new Pythia8::Pythia("../xmldoc",false));
  settings.ConfigurePythiaSignal(pythiaHS.get());
  settings.ConfigurePythiaPileup(pythiaPU.get());
    
  // VBFTaggingAnalysis
  VBFTaggingAnalysis analysis(pythiaHS.get(),pythiaPU.get(), settings);
  analysis.Initialize();
  
  // Event loop
  cout << "Progress:" << endl;
  for (int iev = 0; iev < settings.nEvents; iev++) {
    if (iev%20==0)
      cout << "\tCurrent: " << iev << endl;
    analysis.AnalyzeEvent(iev, settings.pileup);
  }  
  cout << "VBFTagging Analysis Complete!" << endl;
  
  return 0;
}

