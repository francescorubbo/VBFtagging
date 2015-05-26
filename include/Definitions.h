//-*-c++-*-

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define _USE_MATH_DEFINES

#include <utility>
#include <vector>
#include <memory>
#include <random>
#include <set>
#include <math.h>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

const double PI = M_PI;

typedef std::vector<fastjet::PseudoJet> JetVector;

#endif
