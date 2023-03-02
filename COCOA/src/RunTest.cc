#include "RunTest.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "Config_reader_var.hh"

#include <math.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <string>

int RunTest( test_t type ) {

    if ( !type )
	return 0;

    if ( type == kJETS ) {
	G4cout << "Warning. Test type 'Jets' not implemented yet." << G4endl;
	return -1;
    }

    Config_reader_var config_var = Config_reader_var::GetInstance();

    TFile rootFile(config_var.Output_file_path.c_str(), "READ");
    TTree *outTreeLowRes = (TTree*)rootFile.Get("Out_Tree");
    std::vector<float> *cell_energies            = nullptr;
    std::vector<float> *particle_energies        = nullptr;
    TBranch            *branch_cell_energies     = nullptr;
    TBranch            *branch_particle_energies = nullptr;
    outTreeLowRes->SetBranchAddress("cell_e", &cell_energies, &branch_cell_energies);
    outTreeLowRes->SetBranchAddress("particle_e", &particle_energies, &branch_particle_energies);

    std::vector<float> cellESum_over_particleESum_perEvent;
    int nEvents = outTreeLowRes->GetEntries();
    for ( int iEntry = 0; iEntry < nEvents; ++iEntry ) {
	int treeIndex = outTreeLowRes->LoadTree( iEntry );
	branch_cell_energies->GetEntry( treeIndex );
	branch_particle_energies->GetEntry( treeIndex );
	float cell_e_sum     = accumulate( cell_energies->begin(),
					   cell_energies->end(),
					   0.0 );
	float particle_e_sum = accumulate( particle_energies->begin(),
					   particle_energies->end(),
					   0.0 );
	if ( particle_e_sum != 0.0 )
	  cellESum_over_particleESum_perEvent.push_back( cell_e_sum / particle_e_sum );
    }
    float cellESum_over_particleESum_mean = std::accumulate( cellESum_over_particleESum_perEvent.begin(),
							     cellESum_over_particleESum_perEvent.end(),
							     0.0 ) / nEvents;
    float cellESum_over_particleESum_std  = 0.0;
    for ( float e_sum : cellESum_over_particleESum_perEvent )
	cellESum_over_particleESum_std += pow( e_sum - cellESum_over_particleESum_mean, 2 );
    cellESum_over_particleESum_std = sqrt( cellESum_over_particleESum_std / ( (float)nEvents - 1.0 ) );

    int exit_code = 0;

    float mean_expected  = -1.0;
    float std_expected   = -1.0;
    
    float mean_tolerance = 0.05;
    float std_tolerance  = 0.05;

    std::string test_name = "";

    switch ( type ) {
        case kHADRON:
            mean_expected = 0.5;
            std_expected  = 0.23;
            test_name     = "hadron";
            break;
        case kPI_ZERO:
            mean_expected = 0.77;
            std_expected  = 0.16;
            test_name     = "pi^zero";
            break;
        case kJETS:
            break;
        case kNONE:
            break;
    };

    if ( fabs( cellESum_over_particleESum_mean - mean_expected ) > mean_tolerance )
	exit_code = -1;
    if ( fabs( cellESum_over_particleESum_std - std_expected ) > std_tolerance )
	exit_code = -1;

    G4cout << "==============================================================================" << G4endl;
    G4cout << "" << G4endl;
    G4cout << "     Detector Response Test Run for " << test_name << ", " << G4endl;
    G4cout << "     evaluating cell energy sum over truth particle energy sum : " << G4endl;
    G4cout << "" << G4endl;
    G4cout << "    \tObserved             : " << cellESum_over_particleESum_mean << " +/- " << cellESum_over_particleESum_std << G4endl;
    G4cout << "" << G4endl;
    G4cout << "    \tExpected             : " << mean_expected << " +/- " << std_expected << G4endl;
    G4cout << "" << G4endl;
    G4cout << "    \tMean & Std tolerance : " << mean_tolerance << " & " << std_tolerance << G4endl;
    G4cout << "" << G4endl;
    G4cout << "" << G4endl;
    G4cout << "       ===>>> Detector Response Test Run ";
    if ( exit_code )
	G4cout << "Not ";
    G4cout << "Successful!" << G4endl;
    G4cout << "" << G4endl;
    G4cout << "" << G4endl;
    G4cout << "==============================================================================" << G4endl;

    return exit_code;

}
