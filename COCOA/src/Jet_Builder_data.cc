#include "Jet_Builder_data.hh"

Jet_Builder_data::Jet_Builder_data(std::string prefix)
{
    Prefix = prefix;
}

Jet_Builder_data::~Jet_Builder_data()
{;}

void Jet_Builder_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch(TString(Prefix + "_jet_pt"), "vector<float>", &jet_pt);
    outTree->Branch(TString(Prefix + "_jet_eta"), "vector<float>", &jet_eta);
    outTree->Branch(TString(Prefix + "_jet_phi"), "vector<float>", &jet_phi);
    outTree->Branch(TString(Prefix + "_jet_m"), "vector<float>", &jet_m);
    outTree->Branch(TString(Prefix + "_jet_constituents_jetIndex"), "vector<int>", &jet_constituents_jetIndex);
}

void Jet_Builder_data::fill_cell_var()
{
    jet_constituents_jetIndex = std::vector<int>( n_constituents, -1 );
    int size_jets = jets.size();
    for (int ijet = 0; ijet < size_jets; ijet++)
    {
        jet_pt.push_back(jets.at(ijet).pt());
        jet_eta.push_back(jets.at(ijet).eta());
        float phi = jets.at(ijet).phi();
        if (phi < -M_PI) phi += 2*M_PI;
        if (phi > M_PI) phi -= 2*M_PI;
        jet_phi.push_back(phi);
        jet_m.push_back(jets.at(ijet).m());
        int size_constituents = jets.at(ijet).constituents().size();
        std::vector<float> constituents;
        for (int icon = 0; icon < size_constituents; icon++)
        {
	    int constituent_index                        = jets.at(ijet).constituents().at(icon).user_index();
	    jet_constituents_jetIndex[constituent_index] = ijet;
        }
    }
}

void Jet_Builder_data::clear()
{
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_m.clear();
    jet_constituents_jetIndex.clear();
}
