#ifndef __JET_BUILDER_FUNC_H__
#define __JET_BUILDER_FUNC_H__


#include <vector>
#include <stdio.h>
#include <iostream>
#include <map>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "Config_reader_var.hh"

class Jet_Builder_func
{

public:
    void build_jets(std::vector<fastjet::PseudoJet> &input_jets, std::vector<fastjet::PseudoJet> &output_jets, Jet_parameters jet_par);
    void reset();
    fastjet::ClusterSequence *cs;
    fastjet::JetAlgorithm algorithm(std::string algo);
    fastjet::RecombinationScheme recombination_scheme(std::string reco);
};
#endif // __JET_BUILDER_FUNC_H__