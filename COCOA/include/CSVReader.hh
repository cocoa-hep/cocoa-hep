#ifndef __CSVREADER_H__
#define __CSVREADER_H__


#include "G4UIGAG.hh"
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// #include "DataStorage.hh"

class CSVReader
{
    public:
        CSVReader();
        void read_record(TH1F &hist, std::string filename);
        void read_record(TH3F &hist, std::string filename);
        void FillTrackingConfig() ;
        static CSVReader &GetInstance() 
        {
            static CSVReader csv; 
            return csv; 
        };
		TH1F hits_sigma_A0; 
		TH1F hits_sigma_Z0; 
		TH1F hits_sigma_QP; 
		TH1F hits_sigma_Theta; 
		TH1F hits_sigma_Phi; 
		TH1F hist_recon_eff;
		TH3F hist_pflow_mean_template;
		TH3F hist_pflow_std_template;
};

#endif // __CSVREADER_H__
