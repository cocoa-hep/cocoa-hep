#include "CSVReader.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

CSVReader::CSVReader()
{
    FillTrackingConfig();
}

void CSVReader::read_record(TH1F &hist, std::string filename)
{

    // File pointer
    std::fstream fin;
    // Open an existing file
    fin.open(filename, std::ios::in);
    if (!fin)
    {
        G4cout << "Error!!! There is no " << filename << G4endl;
        exit(1);
    }
    std::string line, word;
    std::vector<std::vector<double>> value;
    std::vector<double> valuerow;
    std::getline(fin, line);
    while (std::getline(fin, line))
    {
        valuerow.clear();
        std::stringstream streamstring;
        streamstring << line;
        while (std::getline(streamstring, word, ','))
        {
            // add all the column data
            // of a row to a vector
            valuerow.push_back(std::stod(word));
        }
        // convert string to integer for comparision
        value.push_back(valuerow);
    }

    int nbins = value.size();
    Double_t *xEdges = new Double_t[nbins + 1];
    for (int i = 0; i < nbins; i++)
    {
        xEdges[i] = value.at(i).at(0);
    }
    xEdges[nbins] = 100000000000; //inf
    const char *char_name = filename.c_str();
    hist = TH1F(char_name, "h1", nbins, xEdges);
    for (int i = 1; i <= nbins; i++)
    {
        hist.SetBinContent(i, value.at(i - 1).at(1));
    }
}

void CSVReader::read_record(TH3F &hist, std::string filename)
{

    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(filename, std::ios::in);
    if (!fin)
    {
        G4cout << "Error!!! There is no " << filename << G4endl;
        exit(0);
    }
    std::string line, word;
    std::vector<std::vector<double>> value;
    std::vector<double> valuerow;
    std::getline(fin, line);
    while (std::getline(fin, line))
    {
        valuerow.clear();
        std::stringstream streamstring;
        streamstring << line;
        while (std::getline(streamstring, word, ','))
        {
            // add all the column data
            // of a row to a vector
            valuerow.push_back(std::stod(word));
        }
        // convert string to integer for comparision

        value.push_back(valuerow);
    }

    int nbins = value.size();

    int num_x_bins = 1;
    int num_y_bins = 1;
    int num_z_bins = 1;
    float buff_x = value.at(0).at(0);
    float buff_y = value.at(0).at(1);
    float buff_z = value.at(0).at(2);
    for (int i = 0; i < nbins; i++)
    {
        if (buff_x < value.at(i).at(0))
            num_x_bins++;
        else if (buff_x > value.at(i).at(0))
            num_x_bins = 1;
        if (buff_y < value.at(i).at(1))
            num_y_bins++;
        else if (buff_y > value.at(i).at(1))
            num_y_bins = 1;
        if (buff_z < value.at(i).at(2))
            num_z_bins++;
        else if (buff_z > value.at(i).at(2))
            num_z_bins = 1;
        buff_x = value.at(i).at(0);
        buff_y = value.at(i).at(1);
        buff_z = value.at(i).at(2);
    }
    Double_t *xEdges = new Double_t[num_x_bins + 1];
    Double_t *yEdges = new Double_t[num_y_bins + 1];
    Double_t *zEdges = new Double_t[num_z_bins + 1];
    int i_x = 0;
    int i_y = 0;
    int i_z = 0;
    xEdges[0] = value.at(0).at(0);
    yEdges[0] = value.at(0).at(1);
    zEdges[0] = value.at(0).at(2);
    for (int i = 0; i < nbins; i++)
    {
        if (xEdges[i_x] < value.at(i).at(0))
            i_x++;
        else if (xEdges[i_x] > value.at(i).at(0))
            i_x = 0;
        if (yEdges[i_y] < value.at(i).at(1))
            i_y++;
        else if (yEdges[i_y] > value.at(i).at(1))
            i_y = 0;
        if (zEdges[i_z] < value.at(i).at(2))
            i_z++;
        else if (zEdges[i_z] > value.at(i).at(2))
            i_z = 0;
        xEdges[i_x] = value.at(i).at(0);
        yEdges[i_y] = value.at(i).at(1);
        zEdges[i_z] = value.at(i).at(2);
    }
    xEdges[num_x_bins] = 100000000000;
    yEdges[num_y_bins] = 100000000000;
    zEdges[num_z_bins] = 100000000000;

    const char *char_name = filename.c_str();
    hist = TH3F(char_name, "h3", num_x_bins, xEdges, num_y_bins, yEdges, num_z_bins, zEdges);
    buff_x = value.at(0).at(0);
    buff_y = value.at(0).at(1);
    buff_z = value.at(0).at(2);
    num_x_bins = 1;
    num_y_bins = 1;
    num_z_bins = 1;
    for (int i = 1; i <= nbins; i++)
    {
        if (buff_x < value.at(i - 1).at(0))
            num_x_bins++;
        else if (buff_x > value.at(i - 1).at(0))
            num_x_bins = 1;
        buff_x = value.at(i - 1).at(0);
        if (buff_y < value.at(i - 1).at(1))
            num_y_bins++;
        else if (buff_y > value.at(i - 1).at(1))
            num_y_bins = 1;
        if (buff_z < value.at(i - 1).at(2))
            num_z_bins++;
        else if (buff_z > value.at(i - 1).at(2))
            num_z_bins = 1;
        buff_x = value.at(i - 1).at(0);
        buff_y = value.at(i - 1).at(1);
        buff_z = value.at(i - 1).at(2);
        hist.SetBinContent(num_x_bins, num_y_bins, num_z_bins, value.at(i - 1).at(3));
    }
    delete[] xEdges;
    delete[] yEdges;
    delete[] zEdges;
}

void CSVReader::FillTrackingConfig()
{
    // auto runData
    // 	= static_cast<DataStorage*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    read_record(hits_sigma_A0, "./tracking_configuration/d0_resolution.csv");
    read_record(hits_sigma_Z0, "./tracking_configuration/z0_resolution.csv");
    read_record(hits_sigma_QP, "./tracking_configuration/qp_resolution.csv");
    read_record(hits_sigma_Theta, "./tracking_configuration/theta_resolution.csv");
    read_record(hits_sigma_Phi, "./tracking_configuration/phi_resolution.csv");
    read_record(hist_recon_eff, "./tracking_configuration/reconstruction_efficiency.csv");
    read_record(hist_pflow_mean_template, "./pflow_configuration/pflow_mean_template.csv");
    read_record(hist_pflow_std_template, "./pflow_configuration/pflow_std_template.csv");
}
