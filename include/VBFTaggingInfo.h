#ifndef VBFTAGGINGINFO_H
#define VBFTAGGINGINFO_H

#include "fastjet/PseudoJet.hh"

using namespace fastjet;

class VBFTaggingInfo : public PseudoJet::UserInfoBase{
 public:
 
 VBFTaggingInfo(const int & pdg_id_in,
		const int & pythia_id_in,  
		const double & pv_in, 
		const bool & pileup_in)
   :_pdg_id(pdg_id_in),
    _pythia_id(pythia_id_in), 
    _pv(pv_in),
    _pileup(pileup_in){}
  
  int pdg_id() const { return _pdg_id;}
  int pythia_id() const {return _pythia_id;}
  bool pileup() const { return _pileup;}
  double pv() const { return _pv;}  

 protected:
  int _pdg_id;         // the associated pdg id
  int _pythia_id;  // index in pythia.event
  double _pv;  // the particle pv
  bool _pileup; //true if pileup, false if truth
};

#endif
