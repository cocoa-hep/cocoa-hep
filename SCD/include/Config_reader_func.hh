#ifndef __CONFIG_READER_FUNC_H__
#define __CONFIG_READER_FUNC_H__

// #include "Config_reader_var.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh" // used for G4int, double
#include <string>
#include <json/json.h>
#include "G4NistManager.hh"
#include "Config_reader_var.hh"


class Config_reader_func
{
    public:
        Config_reader_func(std::string path, Config_reader_var &config_var);
    private:
        Json::Value configs;
        void Fill_1D_vector(const Json::Value& Json_list, std::vector<std::vector<long double>>& vec);
        void Fill_1D_vector(const Json::Value &Json_list, std::vector<float> &vec);
        void Fill_1D_vector(const Json::Value &Json_list, std::vector<int> &vec);
        void Fill_1D_vector(const Json::Value& Json_list, std::vector<std::vector<int>>& vec);
        void Fill_1D_vector(const Json::Value& Json_list, std::vector<std::string>& vec);
        void Fill_2D_vector(const Json::Value& Json_list, std::vector< std::vector<long double> >& vec);
        void Fill_2D_vector(const Json::Value& Json_list, std::vector< std::vector<int> >& vec);
        G4Material *Material_build(std::string name);

};






#endif // __CONFIG_READER_FUNC_H__