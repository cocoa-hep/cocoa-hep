#include "ReduceResolution.hh"


ReduceResolution::ReduceResolution(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low)
{
    init_low = 0;
    init_high = 0;
    sum_pixels(CellArray_High, CellArray_Low, config_var.low_resolution.number_of_pixels_ECAL, config_var.high_resolution.number_of_pixels_ECAL);
    sum_pixels(CellArray_High, CellArray_Low, config_var.low_resolution.number_of_pixels_HCAL, config_var.high_resolution.number_of_pixels_HCAL);
}

ReduceResolution::ReduceResolution() 
{;}

ReduceResolution::~ReduceResolution()
{
    //
}


void ReduceResolution::sum_pixels(std::vector<std::vector<std::vector<Cell>>> &CellArray_High,
                                  std::vector<std::vector<std::vector<Cell>>> &CellArray_Low,
                                  std::vector<std::vector<int>> Low_pixel, std::vector<std::vector<int>> High_pixel)
{

    for (int ilow_lay = 0; ilow_lay < (int)Low_pixel.size(); ilow_lay++)
    {
        for (int ilow_eta = 0; ilow_eta < Low_pixel.at(ilow_lay).at(0); ilow_eta++)
        {
            for (int ilow_phi = 0; ilow_phi < Low_pixel.at(ilow_lay).at(0); ilow_phi++)
            {
                // CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).set_noise_signal(gRandom->Gaus(0., config_var.low_resolution.layer_noise_ECAL.at(ilow_lay).at(0)) * MeV);
                std::vector<Particle_dep_in_cell> ParticlesInClust;
                double CellTotSignal = 0; //NoiseArray[it_x][it_y];
                double CellChSignal = 0;  //CellTotSignal;
                double CellNuSignal = 0;  //CellTotSignal;
                std::vector<std::vector<int>> low_cell_link;
                std::vector<std::vector<int>> high_cell_link= {{ilow_lay + init_low, ilow_eta, ilow_phi}};
                for (int ihigh_lay = 0; ihigh_lay < (int)High_pixel.at(ilow_lay).size(); ihigh_lay++)
                {
                    int scale = High_pixel.at(ilow_lay).at(ihigh_lay) / Low_pixel.at(ilow_lay).at(0);
                    for (int ihigh_eta = scale * ilow_eta; ihigh_eta < scale * (ilow_eta + 1); ihigh_eta++)
                    {
                        for (int ihigh_phi = scale * ilow_phi; ihigh_phi < scale * (ilow_phi + 1); ihigh_phi++)
                        {
                            low_cell_link.push_back({init_high + ihigh_lay, ihigh_eta, ihigh_phi});
                            CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).set_cell_link_superres(high_cell_link);

                            CellTotSignal += CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_total_energy();
                            CellChSignal += CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_charge_energy();
                            CellNuSignal += CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_neutral_energy();
                            CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).add_charge_energy(CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_charge_energy());
                            CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).add_neutral_energy(CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_neutral_energy());
                            //
                            std::vector<Particle_dep_in_cell> prtls;
                            CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_particles(prtls);
                            int sizeprtcl = prtls.size();
                            for (int prtclN = 0; prtclN < sizeprtcl; prtclN++)
                            {
                                CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).add_particle(prtls.at(prtclN));
                            }
			    std::vector<Particle_dep_in_cell> conv_el_all;
                            CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).get_conv_electrons( conv_el_all );
			    size_t n_conv_el = conv_el_all.size();
			    for( size_t iConvEl = 0; iConvEl < n_conv_el; ++iConvEl )
				CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).add_particle( conv_el_all[iConvEl], true );
			}
                    }
                }
                CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).set_cell_link_superres(low_cell_link);
            }
        }
        init_high += High_pixel.at(ilow_lay).size();
    }
    init_low = Low_pixel.size();
    apply_noise(CellArray_Low);
}

void ReduceResolution::link_apply(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low)
{
    init_low = 0;
    init_high = 0;
    topo_label_apply(CellArray_High, CellArray_Low, config_var.low_resolution.number_of_pixels_ECAL, config_var.high_resolution.number_of_pixels_ECAL);
    topo_label_apply(CellArray_High, CellArray_Low, config_var.low_resolution.number_of_pixels_HCAL, config_var.high_resolution.number_of_pixels_HCAL);
}

void ReduceResolution::apply_noise(std::vector<std::vector<std::vector<Cell>>> &CellArray) 
{
    for (int ilay = 0; ilay < (int)CellArray.size(); ilay++)
    {
        for (int ieta = 0; ieta < (int)CellArray.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)CellArray.at(ilay).at(ieta).size(); iphi++)
            {
                CellArray.at(ilay).at(ieta).at(iphi).set_noise_signal(gRandom->Gaus(0., config_var.low_resolution.layer_noise.at(ilay)) * MeV);
            }
        }
    }
}

void ReduceResolution::topo_label_apply(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low,
                                        std::vector<std::vector<int>> Low_pixel, std::vector<std::vector<int>> High_pixel)
{
    for (int ilow_lay = 0; ilow_lay < (int)Low_pixel.size(); ilow_lay++)
    {
        for (int ilow_eta = 0; ilow_eta < Low_pixel.at(ilow_lay).at(0); ilow_eta++)
        {
            for (int ilow_phi = 0; ilow_phi < Low_pixel.at(ilow_lay).at(0); ilow_phi++)
            {
                int label = CellArray_Low.at(ilow_lay + init_low).at(ilow_eta).at(ilow_phi).get_label();
                for (int ihigh_lay = 0; ihigh_lay < (int)High_pixel.at(ilow_lay).size(); ihigh_lay++)
                {
                    int scale = High_pixel.at(ilow_lay).at(ihigh_lay) / Low_pixel.at(ilow_lay).at(0);
                    for (int ihigh_eta = scale * ilow_eta; ihigh_eta < scale * (ilow_eta + 1); ihigh_eta++)
                    {
                        for (int ihigh_phi = scale * ilow_phi; ihigh_phi < scale * (ilow_phi + 1); ihigh_phi++)
                        {
                            CellArray_High.at(init_high + ihigh_lay).at(ihigh_eta).at(ihigh_phi).set_label(label);
                        }
                    }
                }
            }
        }
        init_high += High_pixel.at(ilow_lay).size();
    }
    init_low = Low_pixel.size();
}
