#include "Topo_clust_var.hh"

Topo_clust::Topo_clust()
{
    label = 0;
    total_energy = 0;
    EM_energy = 0;
    abs_total_energy = 0;
    noise = 0;
    neutral_energy = 0;
    charge_energy = 0;

    charge = 0;
    eta_com = 0;
    phi_com = 0;
    R_com = 0;
    x_com = 0;
    y_com = 0;
    z_com = 0;
    sigma_eta = 0;
    sigma_phi = 0;
    cos_phi = 0;
    sin_phi = 0;
    n_cell = 0;
    sigma_eta_buf = 0;
    sigma_phi_buf = 0;
    xyz_com.resize(3);
    truth_link = -1;
}
Topo_clust::Topo_clust(int Label, Geometry_definition Geometry)
{
    geometry = Geometry;
    label = Label;
    total_energy = 0;
    EM_energy = 0;
    abs_total_energy = 0;
    noise = 0;
    neutral_energy = 0;
    charge_energy = 0;

    charge = 0;
    eta_com = 0;
    phi_com = 0;
    R_com = 0;
    x_com = 0;
    y_com = 0;
    z_com = 0;
    sigma_eta = 0;
    sigma_phi = 0;
    cos_phi = 0;
    sin_phi = 0;
    n_cell = 0;
    sigma_eta_buf = 0;
    sigma_phi_buf = 0;
    xyz_com.resize(3);
    truth_link = -1;
}

void Topo_clust::add_cell(Cell &cell) //Todo add calculating weight
{
    cell_update(cell, 1, 1);
}
void Topo_clust::subtract_cell(Cell &cell, float fraction)
{
    cell_update(cell, -1, fraction);
}

void Topo_clust::cell_update(Cell &cell ,int add_subtract, float fraction)
{
    n_cell += add_subtract*1;
    if (add_subtract==-1)
        cell_energy = cell.get_total_energy();
    else
        cell_energy = cell.get_pflow_remnant();
    if (cell.get_is_cell_shared())
    {
        if (label == cell.get_1st_local_max_label())
        {
            total_energy += fraction*add_subtract*cell.get_1st_local_max_weight() * cell_energy;
            abs_total_energy += fraction*add_subtract*cell.get_1st_local_max_weight() * fabs(cell_energy);
            if (cell.get_layer() < 3)
                EM_energy += fraction*add_subtract*cell.get_1st_local_max_weight() * cell_energy;
            noise += add_subtract*cell.get_1st_local_max_weight() * cell.get_noise_signal();
            neutral_energy += add_subtract*cell.get_1st_local_max_weight() * cell.get_neutral_energy();
            charge_energy += add_subtract*cell.get_1st_local_max_weight() * cell.get_charge_energy();
            cell.set_truth_label();
            if (particle_type.count((cell.get_truth_label())))
            {
                Energy_Share.at(particle_type.at(cell.get_truth_label())) += add_subtract*cell.get_1st_local_max_weight() * cell.get_total_energy();
            }
            else if (cell.get_truth_label()!=-1)
            {
                Energy_Share.at(4) += add_subtract*cell.get_1st_local_max_weight() * cell.get_total_energy();
            }
            COM_update(cell, fraction*cell.get_1st_local_max_weight(), add_subtract);
        }
        else if (label == cell.get_2nd_local_max_label())
        {
            total_energy += fraction*add_subtract*cell.get_2nd_local_max_weight() * cell_energy;
            abs_total_energy += fraction*add_subtract*cell.get_2nd_local_max_weight() * fabs(cell_energy);
            if (cell.get_layer() < 3)
                EM_energy += fraction*add_subtract*cell.get_2nd_local_max_weight() * cell_energy;
            noise += add_subtract*cell.get_2nd_local_max_weight() * cell.get_noise_signal();
            neutral_energy += add_subtract*cell.get_2nd_local_max_weight() * cell.get_neutral_energy();
            charge_energy += add_subtract*cell.get_2nd_local_max_weight() * cell.get_charge_energy();
            cell.set_truth_label();
            if (particle_type.count((cell.get_truth_label())))
            {
                Energy_Share.at(particle_type.at(cell.get_truth_label())) += add_subtract*cell.get_2nd_local_max_weight() * cell.get_total_energy();
            }
            else if (cell.get_truth_label()!=-1)
            {
                Energy_Share.at(4) += add_subtract*cell.get_2nd_local_max_weight() * cell.get_total_energy();
            }
            COM_update(cell, fraction*cell.get_2nd_local_max_weight(), add_subtract);
        }
        else
        {
            G4cout << "ERORR in Topo_clust::add_cell cell.get_is_cell_shared()=true" << G4endl;
        }
    }
    else
    {
        if (label == cell.get_label())
        {
            total_energy += fraction*add_subtract*cell_energy;
            abs_total_energy += fraction*add_subtract*fabs(cell_energy);
            if (cell.get_layer() < 3)
                EM_energy += fraction*add_subtract*cell_energy;
            noise += add_subtract*cell.get_noise_signal();
            neutral_energy += add_subtract*cell.get_neutral_energy();
            charge_energy += add_subtract*cell.get_charge_energy();
            cell.set_truth_label();
            if (particle_type.count((cell.get_truth_label())))
            {
                Energy_Share.at(particle_type.at(cell.get_truth_label())) += add_subtract*cell.get_total_energy();
            }
            else if (cell.get_truth_label()!=-1)
            {
                Energy_Share.at(4) += add_subtract*cell.get_total_energy();
            }
            COM_update(cell, fraction*1, add_subtract);
        }
        else
        {
            G4cout << "ERORR in Topo_clust::add_cell cell.get_is_cell_shared()=false. Cluster label " << label << " Cell label " << cell.get_label() << G4endl;
        }
    }
    
    float theta_com = 2.0 * TMath::ATan(TMath::Exp(-eta_com));
    x_com = R_com * TMath::Cos(phi_com);
    y_com = R_com * TMath::Sin(phi_com);
    z_com = R_com / TMath::Tan(theta_com);
    xyz_com.at(0) = x_com;
    xyz_com.at(1) = y_com;
    xyz_com.at(2) = z_com;
    sigma_eta = sqrt(sigma_eta_buf / n_cell);
    sigma_phi = sqrt(sigma_phi_buf / n_cell);
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    if (config_var.Type_of_running != "PFlow_debug_E_p_template")
    {
        if (sigma_eta < 0.05)
            sigma_eta = 0.05;
        if (sigma_phi < 0.05)
            sigma_phi = 0.05;
    }
    px = total_energy*sin(theta_com)*cos(phi_com);
    py = total_energy*sin(theta_com)*sin(phi_com);
    pz = total_energy*cos(theta_com);
    std::vector<float>::iterator Energy_Share_result = std::max_element(Energy_Share.begin(), Energy_Share.end());
    truth_link = float(std::distance(Energy_Share.begin(), Energy_Share_result));
}


void Topo_clust::COM_update(Cell &cell, float weight, int add_subtract)
{    
    if (abs_total_energy!=0)
    {
        eta_com += add_subtract*fabs(cell_energy) * weight * (cell.get_eta_pos() - eta_com) / abs_total_energy;
        cos_phi += add_subtract*fabs(cell_energy) * weight * (cos(cell.get_phi_pos()) - cos_phi) / abs_total_energy;
        sin_phi += add_subtract*fabs(cell_energy) * weight * (sin(cell.get_phi_pos()) - sin_phi) / abs_total_energy;
        R_com   += add_subtract*fabs(cell_energy) * weight * (geometry.layer_mid_radius_flatten.at(cell.get_layer()) - R_com)/abs_total_energy;
        phi_com = atan2f(sin_phi, cos_phi);
        float argumEta = pow(cell.get_eta_pos() - eta_com, 2);
        sigma_eta_buf += add_subtract*argumEta;//!  sometimes is negative

        float DistancePhi = fabs(cell.get_phi_pos() - phi_com);
        if (DistancePhi > M_PI)
        {
            DistancePhi -= 2 * M_PI;
        }
        float argumPhi = pow(DistancePhi, 2);
        sigma_phi_buf += add_subtract*argumPhi;//!  sometimes is negative
    }
}
