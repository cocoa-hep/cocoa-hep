#include "Config_reader_var.hh"


Config_reader_var::Config_reader_var()  :

    run_hadron_test( false ),
    run_piZero_test( false ),
    run_jets_test( false )
    
{
    Output_file_path = "PFlowNtupleFile_QCD.root";
    Type_of_running  = "Standard";
    Save_truth_particle_graph = false;
    Use_high_granularity = false;
    use_inner_detector = true;
    r_inn_calo = 1500;
    Layer_gap = 8;
    fieldValue = 0.0;
    Material_ECAL = NULL;
    Material_HCAL = NULL;
}
