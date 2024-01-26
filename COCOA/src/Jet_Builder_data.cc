#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "Jet_Builder_data.hh"
#include "HepMCG4Interface.hh"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenParticle_fwd.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex_fwd.h"
#include "HepMC3/GenVertex.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"

#include <algorithm>
#include <cmath>

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
    if ( Prefix == "true" )
	outTree->Branch(TString(Prefix + "_jet_truth_id"), "vector<int>", &jet_truth_id);
}

void Jet_Builder_data::fill_cell_var( float truth_jet_radius )
{
    
    //
    // Looking for truth particle origin candidates:
    // u, d, s, c, b, e, mu, tau, photon, gluon
    //
    std::vector<int>                         prod_vertex_mother_abs_pids = { 1, 2, 3, 4, 5, 11, 13, 15, 22, 21 };
    std::vector<HepMC3::ConstGenParticlePtr> truth_candidates;
    HepMC3::GenEvent*                        hepmcevt = nullptr;
    if ( Prefix == "true" ) {
	PrimaryGeneratorAction *gen_action = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	G4String                genName    = gen_action->GetGeneratorName();
	if ( genName == "hepmc" || genName == "pythia8" ) {
	    HepMC3::GenEvent* hepmcevt = ((HepMCG4Interface*)gen_action->GetGenerator())->GetHepMCGenEvent();
	    for ( HepMC3::ConstGenParticlePtr prt : hepmcevt->particles() ) {
		if ( std::find( prod_vertex_mother_abs_pids.begin(),
				prod_vertex_mother_abs_pids.end(),
				prt->abs_pid() ) != prod_vertex_mother_abs_pids.end() )
		    truth_candidates.push_back( prt );
	    }
	}
    }
    
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
	if ( Prefix == "true" ) {
	    // Find best match among truth particle candidates: within dR < R_jet and best pT fit
	    // ToDo: "overlap removal" - preference for some particles, e.g. e before tau, photon before gluon.
	    float    dR               = 0.0;
	    float    pT_abs_rel_best  = 0.0;
	    int      abs_pid_best_fit = 0; // default pid 0 : no matching truth particle found
	    TVector3 jet_p3( jets[ijet].px(), jets[ijet].py(), jets[ijet].pz() );
	    for ( HepMC3::ConstGenParticlePtr prt : truth_candidates ) {
		TVector3 particle_p3( prt->momentum().px(),  prt->momentum().py(),  prt->momentum().pz() );
		dR = ROOT::Math::VectorUtil::DeltaR( jet_p3, particle_p3 );
		if ( dR > truth_jet_radius )
		    continue;
		float pT_abs_rel = abs( jet_p3.Pt() / particle_p3.Pt() - 1.0 );
		if ( abs_pid_best_fit == 0 || pT_abs_rel < pT_abs_rel_best ) {
		    pT_abs_rel_best  = pT_abs_rel;
		    abs_pid_best_fit = prt->abs_pid();
		}
	    }
	    jet_truth_id.push_back( abs_pid_best_fit );	    
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
    jet_truth_id.clear();
}
