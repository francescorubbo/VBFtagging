#include "Configuration.h"

Configuration::Configuration(int argc, char* argv[]){

    po::options_description gen_desc("Allowed options");
    gen_desc.add_options()
      ("help", "produce help message")
      ("Debug",     po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile",   po::value<string>(&outName)->default_value("VBFTagging.root"), "output file name")
      ("Seed",      po::value<int>(&seed)->default_value(-1), "Seed. -1 means random seed");

    po::options_description sim_desc("Simulation Settings");
    sim_desc.add_options()
      ("NEvents",   po::value<int>(&nEvents)->default_value(1) ,    "Number of Events ")
      ("Pileup",    po::value<int>(&pileup)->default_value(80), "Number of Additional Interactions")
      ("Proc",      po::value<int>(&proc)->default_value(4), "Process:\n - 1: Z'T->ttbar\n - 2: W'->WZ+lept\n - 3: W'->WZ+had\n - 4: QCD")
      ("pThatMin",  po::value<float>(&pThatmin)->default_value(100), "pThatMin for QCD")
      ("pThatMax",  po::value<float>(&pThatmax)->default_value(500), "pThatMax for QCD")
      ("BosonMass", po::value<float>(&boson_mass)->default_value(1500), "Z' or W' mass in GeV");
    
    po::options_description desc;
    desc.add(gen_desc).add(sim_desc);
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
	exit(0);
    }

    print();

    seed=getSeed(seed);    
}

void Configuration::print(){

  cout << endl << "=================================================================" << endl;
  cout << "=                        VBFTagging Analysis                        =" << endl;
  cout << "=================================================================" << endl << endl;

  cout << "Settings:" << endl;
  for (po::variables_map::const_iterator itr=vm.begin();itr != vm.end();++itr){
    printf("%15s\t",itr->first.c_str());
    
    try { 
      cout << "= " << itr->second.as<double>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<float>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<int>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<std::string>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
  }
}

void Configuration::ConfigurePythiaSignal(Pythia8::Pythia* hs){

    hs->readString("Print:quiet=on");
    hs->readString("Random:setSeed = on"); 
    std::stringstream ss; 
    ss << "Random:seed = " << seed;
    hs->readString(ss.str());

   if(proc ==1){
     std::stringstream bosonmass_str; 
     bosonmass_str<< "32:m0=" << boson_mass;
     hs->readString(bosonmass_str.str());
     hs->readString("NewGaugeBoson:ffbar2gmZZprime= on");
     hs->readString("Zprime:gmZmode=3");
     hs->readString("32:onMode = off");
     hs->readString("32:onIfAny = 6");
     hs->readString("24:onMode = off");
     hs->readString("24:onIfAny = 1 2 3 4");
     hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line! 
   }else if(proc ==2){
      std::stringstream bosonmass_str; 
      bosonmass_str<< "34:m0=" << boson_mass;
      
      hs->readString(bosonmass_str.str());
      hs->readString("NewGaugeBoson:ffbar2Wprime = on");
      hs->readString("Wprime:coup2WZ=1");
      hs->readString("34:onMode = off");
      hs->readString("34:onIfAny = 23 24");
      hs->readString("24:onMode = off");
      hs->readString("24:onIfAny = 1 2 3 4");
      hs->readString("23:onMode = off");
      hs->readString("23:onIfAny = 12");
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 3){
      std::stringstream bosonmass_str; 
      bosonmass_str<< "34:m0=" << boson_mass;
      hs->readString(bosonmass_str.str());
      hs->readString("NewGaugeBoson:ffbar2Wprime = on");
      hs->readString("Wprime:coup2WZ=1");
      hs->readString("34:onMode = off");
      hs->readString("34:onIfAny = 23 24");
      hs->readString("24:onMode = off");
      hs->readString("24:onIfAny = 11 12");
      hs->readString("23:onMode = off");
      hs->readString("23:onIfAny = 1 2 3 4 5");
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 4){ 
      hs->readString("HardQCD:all = on");
      std::stringstream ptHatMin;
      std::stringstream ptHatMax;
      ptHatMin << "PhaseSpace:pTHatMin  =" << pThatmin;
      ptHatMax << "PhaseSpace:pTHatMax  =" << pThatmax;
      hs->readString(ptHatMin.str());
      hs->readString(ptHatMax.str());
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }
   else{ 
     throw std::invalid_argument("received invalid 'process'");
   }
   return;
}

void Configuration::ConfigurePythiaPileup(Pythia8::Pythia* pu){
   //Setup the pileup
   pu->readString("Random:setSeed = on");   
   std::stringstream ss;
   ss.clear(); 
   ss.str(""); 
   ss << "Random:seed = " << seed+1; 
   pu->readString(ss.str());
   pu->readString("Print:quiet=on");
   pu->readString("SoftQCD:nonDiffractive = on");
   pu->readString("HardQCD:all = off");
   pu->readString("PhaseSpace:pTHatMin  = .1");
   pu->readString("PhaseSpace:pTHatMax  = 20000");
   pu->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */);
   return;
}

int Configuration::getSeed(int seed){                                                      
  if (seed > -1) return seed;
  int timeSeed = time(NULL);                                                                 
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}
