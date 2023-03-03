#include "Topo_clust_func.hh"

int modulo(int value, int mo)
{
    int mod = value % (int)mo;
    if (value < 0)
    {
        mod += mo;
    }
    return mod;
}

Topo_clust_func::Topo_clust_func(std::vector<std::vector<std::vector<Cell>>> &cells_array, Geometry_definition Geometry, Topological_clustering_var Topo_Config ,std::string class_of_debug) : Cells_Array(cells_array)
{
    topo_config = Topo_Config;
    geometry = Geometry;
    Nlayers = Cells_Array.size();
    int n = 0;
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                n += 1;
                Cell &local_cell = Cells_Array.at(ilay).at(ieta).at(iphi);
                if (class_of_debug == "neutral")
                {
                    local_cell.add_charge_energy(-1 * local_cell.get_charge_energy());
                }
                else if (class_of_debug == "charge")
                {
                    local_cell.add_neutral_energy(-1 * local_cell.get_neutral_energy());
                }
                else if (class_of_debug == "noise")
                {
                    local_cell.add_charge_energy(-1 * local_cell.get_charge_energy());
                    local_cell.add_neutral_energy(-1 * local_cell.get_neutral_energy());
                }

                if(topo_config.cluster_negative_energy_cells)
                {
                    local_cell.set_signal_to_noise(fabs(local_cell.get_total_energy() / geometry.layer_noise.at(ilay))); //sigma not signed
                }
                else
                {
                    local_cell.set_signal_to_noise(local_cell.get_total_energy() / geometry.layer_noise.at(ilay)); //sigma signed (negative energy cells will not be 4-2-0 clustered)
                }
            }
        }
    }
}

void Topo_clust_func::topoclustering(std::vector<Topo_clust> &topo_clusts_list)
{
    std::vector<Cell *> seed_list;
    std::vector<Cell *> local_max_list;
    std::vector<Cell *> share_list;
    seed_list.reserve(50000000);
    local_max_list.reserve(50000000);
    share_list.reserve(50000000);

    find_seed_cells(seed_list);
    cluster_maker(seed_list);
    find_local_max(local_max_list);
    cluster_split(local_max_list, share_list);
    fill_clusters_list(share_list, topo_clusts_list);
}

void Topo_clust_func::insert_to_order_cell_list(Cell &cell, std::vector<Cell *> &cells_vector)
{
    int size_cells_vector = cells_vector.size();
    if (size_cells_vector == 0)
    {
        cell.set_label(1);
        cells_vector.push_back(&cell);
    }
    else
    {
        Cell *last_seed = cells_vector.back();
        if (last_seed->get_signal_to_noise() > cell.get_signal_to_noise())
        {
            cell.set_label(last_seed->get_label() + 1);
            cells_vector.push_back(&cell);
        }
        else
        {
            for (int iseed = size_cells_vector - 1; iseed >= 0; iseed--)
            {
                int label = cells_vector.at(iseed)->get_label();
                if (cells_vector.at(iseed)->get_signal_to_noise() >= cell.get_signal_to_noise())
                {
                    cell.set_label(label + 1);
                    cells_vector.emplace(cells_vector.begin() + iseed + 1, &cell);
                    break;
                }
                else
                    cells_vector.at(iseed)->set_label(label + 1);
                if (iseed == 0)
                {
                    cell.set_label(1);
                    cells_vector.emplace(cells_vector.begin(), &cell);
                }
            }
        }
    }
}

void Topo_clust_func::find_seed_cells(std::vector<Cell *> &seeds)
{
    int n = 0;
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {		
                Cell &local_cell = Cells_Array.at(ilay).at(ieta).at(iphi);
                float local_sigma = local_cell.get_signal_to_noise();
                if (local_sigma > topo_config.sigma_threshold_for_seed_cells)
                {
                    n += 1;
                    insert_to_order_cell_list(local_cell, seeds);
                }
            }
        }
    }
}

int Topo_clust_func::number_of_labels()
{
    int num_labels = 0;
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                int label = Cells_Array.at(ilay).at(ieta).at(iphi).get_label();
                if (label > num_labels)
                    num_labels = label;
            }
        }
    }
    return num_labels;
}

void Topo_clust_func::cluster_maker(std::vector<Cell *> &seeds)
{
    int size_seeds = seeds.size();

    int n = 0;
    number_of_labels();
    while (seeds.size() > 0)
    {
        n++;
        cluster_maker_neighbor(seeds, size_seeds);
        size_seeds--;
        if (size_seeds == 0)
            size_seeds = seeds.size();
    }
    int num_of_labels = number_of_labels();

    //* Deleting Topoclusters with energy below zero (net negative energy)
    std::vector<float> clusters_energy(num_of_labels, 0.0);
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                if (Cells_Array.at(ilay).at(ieta).at(iphi).get_label() != 0)
                {
                    clusters_energy.at(Cells_Array.at(ilay).at(ieta).at(iphi).get_label() - 1) += Cells_Array.at(ilay).at(ieta).at(iphi).get_total_energy();
                }
            }
        }
    }
    for (int itopo = 0; itopo < num_of_labels; itopo++)
    {
        if (clusters_energy.at(itopo) <= 0)
        {
            for (int ilay = 0; ilay < Nlayers; ilay++)
            {
                for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
                {
                    for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
                    {
                        int label = Cells_Array.at(ilay).at(ieta).at(iphi).get_label();
                        if (label == itopo + 1 || label == 0)
                            Cells_Array.at(ilay).at(ieta).at(iphi).set_label(0);
                        else if (label > itopo + 1)
                            Cells_Array.at(ilay).at(ieta).at(iphi).set_label(Cells_Array.at(ilay).at(ieta).at(iphi).get_label() - 1);
                    }
                }
            }
            clusters_energy.erase(clusters_energy.begin() + itopo);
            num_of_labels -= 1;
            itopo -= 1;
        }
    }
}

void Topo_clust_func::cluster_maker_neighbor(std::vector<Cell *> &seeds, int list_size)
{
    Cell *cell = seeds.front();
    std::vector<float> layer_step(Nlayers - 1, 1);
    for (int ilay = 0; ilay < Nlayers - 1; ilay++)
    {
        layer_step.at(ilay) = ( (float)Cells_Array.at(ilay).size() / (float)Cells_Array.at(ilay + 1).size());
    }
    cell->set_size_of_seeds_list(list_size);

    for (int ilay_step = -1; ilay_step < 2; ilay_step += 2)
    {
        int next_layer = cell->get_layer() + ilay_step;
        if (next_layer >= 0 && next_layer < Nlayers)
        {
            if (ilay_step == 1) //Todo: try if next layer has high gran. by if (layer_step.at(next_layer-1] > 1)
            {
                Cell &neighbor_cell = Cells_Array.at(next_layer).at((int)floor(((float)cell->get_eta()) / layer_step.at(next_layer - 1))).at((int)floor(((float)cell->get_phi()) / layer_step.at(next_layer - 1)));
                cluster_maker_add_cell(*cell, neighbor_cell, seeds);
            }
            else //* ilay_step==-1
            {
                for (int layer_stepStepEta = 0; layer_stepStepEta < layer_step.at(next_layer); layer_stepStepEta++)
                {
                    for (int layer_stepStepPhi = 0; layer_stepStepPhi < layer_step.at(next_layer); layer_stepStepPhi++)
                    {
                        Cell &neighbor_cell = Cells_Array.at(next_layer).at(cell->get_eta() * layer_step.at(next_layer) + layer_stepStepEta).at(cell->get_phi() * layer_step.at(next_layer) + layer_stepStepPhi);
                        cluster_maker_add_cell(*cell, neighbor_cell, seeds);
                    }
                }
            }
        }
    }

    for (int deta = -1; deta < 2; deta++)
    {
        for (int dphi = -1; dphi < 2; dphi++)
        {
            int step_eta = cell->get_eta() + deta;
            int step_phi = cell->get_phi() + dphi;
            if (step_eta >= 0 && step_eta < (int)Cells_Array.at(cell->get_layer()).size() && (deta != 0 || dphi != 0))
            {
                Cell &neighbor_cell = Cells_Array.at(cell->get_layer()).at(step_eta).at(modulo(step_phi, Cells_Array.at(cell->get_layer()).at(step_eta).size()));
                cluster_maker_add_cell(*cell, neighbor_cell, seeds);
            }
        }
    }
    seeds.erase(seeds.begin());
}

void Topo_clust_func::cluster_maker_add_cell(Cell &cell, Cell &neighbor_cell, std::vector<Cell *> &seeds)
{

    int local_list_size_4sigma = cell.get_size_of_seeds_list();
    float neighbor_sigma = neighbor_cell.get_signal_to_noise();
    int neighbor_label = neighbor_cell.get_label();

    if (neighbor_sigma > topo_config.sigma_threshold_for_neighboring_cells) //Todo
    {

        if (neighbor_label > cell.get_label())
        {

            global_label_change(neighbor_label, cell.get_label());
        }
        else if (neighbor_label != 0 && neighbor_label < cell.get_label())
        {

            global_label_change(cell.get_label(), neighbor_label);
        }
        else if (neighbor_label == 0)
        {
            neighbor_cell.set_label(cell.get_label());
            int size_seeds = seeds.size();

            if (seeds.back()->get_signal_to_noise() > neighbor_sigma)
            {
                seeds.push_back(&neighbor_cell);
            }
            else
            {
                if (size_seeds == cell.get_size_of_seeds_list())
                {
                    seeds.push_back(&neighbor_cell);
                }
                else
                {
                    for (int iseed = size_seeds - 1; iseed >= local_list_size_4sigma; iseed--)
                    {
                        if (seeds.at(iseed)->get_signal_to_noise() > neighbor_sigma)
                        {
                            seeds.emplace(seeds.begin() + iseed + 1, &neighbor_cell);
                            break;
                        }
                        if (iseed == local_list_size_4sigma)
                        {
                            seeds.insert(seeds.begin() + iseed, &neighbor_cell);
                        }
                    }
                }
            }
        }
    }
    else if (neighbor_sigma > topo_config.sigma_threshold_for_last_cells)
    {
        if (neighbor_label == 0)
        {
            neighbor_cell.set_label(cell.get_label());
            // Cells_Array.at(neighbor_layer).at(neighbor_eta).at(neighbor_phi).set_label(cell.get_label());
        }
    }
}

void Topo_clust_func::global_label_change(int old_label, int new_label)
{
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                Cell &local_cell = Cells_Array.at(ilay).at(ieta).at(iphi);
                int local_label = local_cell.get_label();
                if (local_label == old_label)
                    local_cell.set_label(new_label);
                else if (local_label > old_label)
                    local_cell.set_label(local_label - 1);
            }
        }
    }
}

void Topo_clust_func::find_local_max(std::vector<Cell *> &local_maxs)
{
    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                Cell &cell = Cells_Array.at(ilay).at(ieta).at(iphi);
                float cell_ener = cell.get_total_energy();
                if (is_cell_local_max(cell))
                {
                    int size_local_maxs = local_maxs.size();
                    if (size_local_maxs == 0)
                        local_maxs.push_back(&cell);
                    else
                    {
                        bool check = true;
                        int size_Low_layer_deta_ECAL = geometry.layer_deta_ECAL.size();
                        if (ilay >= size_Low_layer_deta_ECAL)
                        {
                            for (int imax = 0; imax < size_local_maxs; imax++)
                            {
                                if ((((local_maxs.at(imax)->get_eta_pos() >= cell.get_abs_min_eta()) &&
                                      (local_maxs.at(imax)->get_eta_pos() <= cell.get_abs_max_eta())) &&
                                     ((local_maxs.at(imax)->get_phi_pos() >= cell.get_abs_min_phi()) &&
                                      (local_maxs.at(imax)->get_phi_pos() <= cell.get_abs_max_phi()))))
                                {
                                    check = false;
                                }
                            }
                        }
                        if (check)
                        {
                            if (local_maxs.back()->get_total_energy() >= cell_ener)
                                local_maxs.push_back(&cell);
                            else
                            {
                                for (int ilocal = size_local_maxs - 1; ilocal >= 0; ilocal--)
                                {
                                    if (local_maxs.at(ilocal)->get_total_energy() >= cell_ener)
                                    {
                                        local_maxs.emplace(local_maxs.begin() + ilocal + 1, &cell);
                                        break;
                                    }
                                    if (ilocal == 0)
                                    {
                                        local_maxs.emplace(local_maxs.begin(), &cell);
                                    }
                                }
                                
                            }
                        }
                    }
                }
            }
        }
    }
}

bool Topo_clust_func::is_cell_local_max(Cell &cell)
{
    float cell_ener = cell.get_total_energy();
    int num_neighbors = 0;
    bool is_max = true;
    if (cell.get_label()!=0 && cell_ener > topo_config.local_max_seed_energy)
    {
        std::vector<float> layer_step(Nlayers - 1, 1);
        for (int ilay = 0; ilay < Nlayers - 1; ilay++)
        {
            layer_step.at(ilay) = ((float)Cells_Array.at(ilay).size() / (float)Cells_Array.at(ilay + 1).size());
        }

        for (int i = -1; i < 2; i += 2)
        {
            int next_layer = cell.get_layer() + i;
            if (next_layer >= 0 && next_layer < 6)
            {
                if (i == 1)
                {
                    Cell &neighbor_cell = Cells_Array.at(next_layer).at((int)floor(((float)cell.get_eta()) / (layer_step.at(next_layer - 1)))).at((int)floor(cell.get_phi() / (layer_step.at(next_layer - 1))));
                    if (neighbor_cell.get_label() > 0)
                        num_neighbors++;
                    is_max = ((neighbor_cell.get_total_energy() < cell_ener) && is_max);
                }
                else
                {
                    for (int layer_stepStepEta = 0; layer_stepStepEta < layer_step.at(next_layer); layer_stepStepEta++)
                    {
                        for (int layer_stepStepPhi = 0; layer_stepStepPhi < layer_step.at(next_layer); layer_stepStepPhi++)
                        {
                            Cell &neighbor_cell = Cells_Array.at(next_layer).at((int)(cell.get_eta() * (layer_step.at(next_layer)) + layer_stepStepEta)).at((int)(cell.get_phi() * (layer_step.at(next_layer)) + layer_stepStepPhi));
                            if (neighbor_cell.get_label() > 0)
                                num_neighbors++;
                            is_max = ((neighbor_cell.get_total_energy() < cell_ener) && is_max);
                        }
                    }
                }
            }
        }
        for (int i = -1; i < 2; i++)
        {
            for (int j = -1; j < 2; j++)
            {
                int step_eta = cell.get_eta() + i;
                int step_phi = cell.get_phi() + j;
                if (step_eta >= 0 && step_eta < (int)Cells_Array.at(cell.get_layer()).size() && (i != 0 || j != 0))
                {
                    Cell &neighbor_cell = Cells_Array.at(cell.get_layer()).at(step_eta).at(modulo(step_phi, Cells_Array.at(cell.get_layer()).size()));
                    if (neighbor_cell.get_label() > 0)
                        num_neighbors++;
                    is_max = ((neighbor_cell.get_total_energy() < cell_ener) && is_max);
                }
            }
        }
    }
    bool output = false;
    if (num_neighbors >= 4 && is_max)
    {
        output = true;
    }
    return output;
}

void Topo_clust_func::cluster_split(std::vector<Cell *> &local_maxs, std::vector<Cell *> &shares)
{
    // * clust_label -- cluster number (instead of i), since we split cluster total number of clusters change
    int clust_label = 0;
    // * SplitClustCount -- amount of new clusters after one iteration of spliting
    int num_of_split_clust = 0;
    int num_share_cells_in_previus_clsters = 0;
    int num_of_labels = number_of_labels();
    for (int itopo = 1; itopo <= num_of_labels; itopo++)
    {
        clust_label += (num_of_split_clust + 1);
        num_of_split_clust = 0;
        std::vector<Cell *> local_max_cells_in_clust;
        int size_local_maxs = local_maxs.size();
        int local_max_label = 0;
        for (int imax = 0; imax < size_local_maxs; imax++)
        {
            if (local_maxs.at(imax)->get_label() == clust_label)
            {
                local_maxs.at(imax)->set_1st_local_max_label(clust_label + local_max_label);
                local_maxs.at(imax)->first_local_max = local_maxs.at(imax);
                local_max_cells_in_clust.push_back(local_maxs.at(imax));
                local_max_label++;
            }
        }
        int size_local_max_cells_in_clust = local_max_cells_in_clust.size();
        if (size_local_max_cells_in_clust > 1)
        {
            num_of_split_clust = local_max_label - 1;

            while (local_max_cells_in_clust.size() > 0)
            {
                cluster_split_neighbor(local_max_cells_in_clust, shares, clust_label, size_local_max_cells_in_clust);
                size_local_max_cells_in_clust--;
                if (size_local_max_cells_in_clust == 0)
                {
                    size_local_max_cells_in_clust = local_max_cells_in_clust.size();
                }
            }
            for (int ishar = num_share_cells_in_previus_clsters; ishar < (int)shares.size(); ishar++)
            {
                cluster_share_neighbor(*shares.at(ishar), shares, clust_label);
            }
            for (int ilay = 0; ilay < Nlayers; ilay++)
            {
                for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
                {
                    for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
                    {
                        int label = Cells_Array.at(ilay).at(ieta).at(iphi).get_label();
                        if (label > clust_label)
                        {
                            Cells_Array.at(ilay).at(ieta).at(iphi).set_label(label + num_of_split_clust);
                        }
                        else if (clust_label == label)
                        {
                            Cells_Array.at(ilay).at(ieta).at(iphi).set_label(Cells_Array.at(ilay).at(ieta).at(iphi).get_1st_local_max_label());
                        }
                    }
                }
            }
            num_share_cells_in_previus_clsters = shares.size();
        }
    }
}

void Topo_clust_func::cluster_share_neighbor(Cell &cell, std::vector<Cell *> &shares, int clust_label)
{
    std::vector<float> layer_step(Nlayers - 1, 1);
    for (int ilay = 0; ilay < Nlayers - 1; ilay++)
    {
        layer_step.at(ilay) = ((float)Cells_Array.at(ilay).size() / (float)Cells_Array.at(ilay + 1).size());
    }
    for (int ilay_step = -1; ilay_step < 2; ilay_step += 2)
    {
        int next_layer = cell.get_layer() + ilay_step;
        if (next_layer >= 0 && next_layer < Nlayers)
        {
            if (ilay_step == 1) //Todo
            {
                int eta_idx = (int)floor(((float)cell.get_eta()) / layer_step.at(next_layer - 1));
                int phi_idx = (int)floor(((float)cell.get_phi()) / layer_step.at(next_layer - 1));
                Cell &neighbor_cell = Cells_Array.at(next_layer).at(eta_idx).at(phi_idx);
                if (neighbor_cell.get_label() == clust_label && neighbor_cell.get_1st_local_max_label() == 0)
                {
                    neighbor_cell.set_is_cell_shared(true);
                    neighbor_cell.set_1st_local_max_label(cell.get_1st_local_max_label());
                    neighbor_cell.set_2nd_local_max_label(cell.get_2nd_local_max_label());
                    neighbor_cell.first_local_max = cell.first_local_max;
                    neighbor_cell.second_local_max = cell.second_local_max;
                    neighbor_cell.set_label(0);
                    shares.push_back(&neighbor_cell);
                }
            }
            else
            {
                for (int layer_step_step_eta = 0; layer_step_step_eta < layer_step.at(next_layer); layer_step_step_eta++)
                {
                    int eta_idx = (int)(cell.get_eta() * layer_step.at(next_layer) + layer_step_step_eta);
                    for (int layer_step_step_phi = 0; layer_step_step_phi < layer_step.at(next_layer); layer_step_step_phi++)
                    {
                        int phi_idx = (int)(cell.get_phi() * layer_step.at(next_layer) + layer_step_step_phi);
                        Cell &neighbor_cell = Cells_Array.at(next_layer).at(eta_idx).at(phi_idx);
                        if (neighbor_cell.get_label() == clust_label && neighbor_cell.get_1st_local_max_label() == 0)
                        {
                            neighbor_cell.set_is_cell_shared(true);
                            neighbor_cell.set_1st_local_max_label(cell.get_1st_local_max_label());
                            neighbor_cell.set_2nd_local_max_label(cell.get_2nd_local_max_label());
                            neighbor_cell.first_local_max = cell.first_local_max;
                            neighbor_cell.second_local_max = cell.second_local_max;
                            neighbor_cell.set_label(0);
                            shares.push_back(&neighbor_cell);
                        }
                    }
                }
            }
        }
        for (int i = -1; i < 2; i++)
        {
            for (int j = -1; j < 2; j++)
            {
                int step_eta = cell.get_eta() + i;
                int step_phi = cell.get_phi() + j;

                if (step_eta >= 0 && step_eta < (int)Cells_Array.at(cell.get_layer()).size() && (i != 0 || j != 0))
                {
                    Cell &neighbor_cell = Cells_Array.at(cell.get_layer()).at(step_eta).at(modulo(step_phi, Cells_Array.at(cell.get_layer()).size()));
                    if (neighbor_cell.get_label() == clust_label && neighbor_cell.get_1st_local_max_label() == 0)
                    {
                        neighbor_cell.set_is_cell_shared(true);
                        neighbor_cell.set_1st_local_max_label(cell.get_1st_local_max_label());
                        neighbor_cell.set_2nd_local_max_label(cell.get_2nd_local_max_label());
                        neighbor_cell.first_local_max = cell.first_local_max;
                        neighbor_cell.second_local_max = cell.second_local_max;
                        neighbor_cell.set_label(0);
                        shares.push_back(&neighbor_cell);
                    }
                }
            }
        }
    }
}

void Topo_clust_func::cluster_split_neighbor(std::vector<Cell *> &local_max_cells_in_clust, std::vector<Cell *> &shares, int clust_label, int list_size)
{
    std::vector<float> layer_step(Nlayers - 1, 1);
    for (int ilay = 0; ilay < Nlayers - 1; ilay++)
    {
        layer_step.at(ilay) = ((float)Cells_Array.at(ilay).size() / (float)Cells_Array.at(ilay + 1).size());
    }
    local_max_cells_in_clust.front()->set_size_of_seeds_list(list_size);
    for (int ilay_step = -1; ilay_step < 2; ilay_step += 2)
    {
        int next_layer = local_max_cells_in_clust.front()->get_layer() + ilay_step;
        if (next_layer >= 0 && next_layer < Nlayers)
        {
            if (ilay_step == 1) //Todo
            {
                int eta_idx = (int)floor(((float)local_max_cells_in_clust.front()->get_eta()) / layer_step.at(next_layer - 1));
                int phi_idx = (int)floor(((float)local_max_cells_in_clust.front()->get_phi()) / layer_step.at(next_layer - 1));
                Cell &neighbor_cell = Cells_Array.at(next_layer).at(eta_idx).at(phi_idx);
                if (neighbor_cell.get_label() == clust_label)
                    cluster_split_add_cell(*local_max_cells_in_clust.front(), local_max_cells_in_clust, shares, neighbor_cell);
            }
            else
            {
                for (int layer_step_step_eta = 0; layer_step_step_eta < layer_step.at(next_layer); layer_step_step_eta++)
                {
                    int eta_idx = (int)(local_max_cells_in_clust.front()->get_eta() * layer_step.at(next_layer) + layer_step_step_eta);
                    for (int layer_step_step_phi = 0; layer_step_step_phi < layer_step.at(next_layer); layer_step_step_phi++)
                    {
                        int phi_idx = (int)(local_max_cells_in_clust.front()->get_phi() * layer_step.at(next_layer) + layer_step_step_phi);
                        Cell &neighbor_cell = Cells_Array.at(next_layer).at(eta_idx).at(phi_idx);
                        if (neighbor_cell.get_label() == clust_label)
                            cluster_split_add_cell(*local_max_cells_in_clust.front(), local_max_cells_in_clust, shares, neighbor_cell);
                    }
                }
            }
        }
    }
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            int step_eta = local_max_cells_in_clust.front()->get_eta() + i;
            int step_phi = local_max_cells_in_clust.front()->get_phi() + j;

            if (step_eta >= 0 && step_eta < (int)Cells_Array.at(local_max_cells_in_clust.front()->get_layer()).size() && (i != 0 || j != 0))
            {
                Cell &neighbor_cell = Cells_Array.at(local_max_cells_in_clust.front()->get_layer()).at(step_eta).at(modulo(step_phi, Cells_Array.at(local_max_cells_in_clust.front()->get_layer()).size()));
                if (neighbor_cell.get_label() == clust_label)
                    cluster_split_add_cell(*local_max_cells_in_clust.front(), local_max_cells_in_clust, shares, neighbor_cell);
            }
        }
    }

    local_max_cells_in_clust.erase(local_max_cells_in_clust.begin());
}

void Topo_clust_func::cluster_split_add_cell(Cell &cell, std::vector<Cell *> &local_max_cells_in_clust, std::vector<Cell *> &shares, Cell &neighbor_cell)
{

    float neighbor_ener = neighbor_cell.get_total_energy();
    float first_clust_label = neighbor_cell.get_1st_local_max_label();
    float second_clust_label = neighbor_cell.get_2nd_local_max_label();
    if (first_clust_label == 0)
    {
        neighbor_cell.first_local_max = cell.first_local_max;
        neighbor_cell.set_1st_local_max_label(cell.get_1st_local_max_label());
        int size_local_maxs = local_max_cells_in_clust.size();
        if (cell.get_total_energy() > neighbor_ener)
        {
            local_max_cells_in_clust.push_back(&neighbor_cell);
        }
        else
        {
            if (size_local_maxs == cell.get_size_of_seeds_list())
                local_max_cells_in_clust.push_back(&neighbor_cell);
            else
            {
                for (int itopo = size_local_maxs - 1; itopo >= cell.get_size_of_seeds_list(); itopo--)
                {
                    if (local_max_cells_in_clust.at(itopo)->get_total_energy() > neighbor_ener)
                    {
                        local_max_cells_in_clust.insert(local_max_cells_in_clust.begin() + itopo + 1, &neighbor_cell);
                        break;
                    }
                    else if (itopo == cell.get_size_of_seeds_list())
                        local_max_cells_in_clust.insert(local_max_cells_in_clust.begin() + itopo, &neighbor_cell);
                }
            }
        }
    }
    else if (first_clust_label != cell.get_1st_local_max_label())
    {
        if (second_clust_label == 0)
        {
            neighbor_cell.second_local_max = cell.first_local_max;
            neighbor_cell.set_2nd_local_max_label(cell.get_1st_local_max_label());
            neighbor_cell.set_is_cell_shared(true);
            int size_local_maxs = local_max_cells_in_clust.size();
            for (int icell = 0; icell < size_local_maxs; icell++)
            {
                if ((local_max_cells_in_clust.at(icell)->get_layer() == neighbor_cell.get_layer()) &&
                    (local_max_cells_in_clust.at(icell)->get_eta() == neighbor_cell.get_eta()) &&
                    (local_max_cells_in_clust.at(icell)->get_phi() == neighbor_cell.get_phi()))
                {
                    local_max_cells_in_clust.erase(local_max_cells_in_clust.begin() + icell);
                    break;
                }
            }
            neighbor_cell.set_label(0);
            shares.push_back(&neighbor_cell);
        }
    }
}

void Topo_clust_func::fill_clusters_list(std::vector<Cell *> &shares, std::vector<Topo_clust> &topo_clusts_list)
{
    float energy_ = 0;
    int num_of_labels = number_of_labels();

    for (int itopo = 0; itopo < num_of_labels; itopo++)
    {
        Topo_clust topoclust(itopo + 1, geometry);
        topo_clusts_list.push_back(topoclust);
    }

    for (int ilay = 0; ilay < Nlayers; ilay++)
    {
        for (int ieta = 0; ieta < (int)Cells_Array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)Cells_Array.at(ilay).at(ieta).size(); iphi++)
            {
                Cell &local_cell = Cells_Array.at(ilay).at(ieta).at(iphi);
                if (local_cell.get_label() != 0)
                {
                    energy_ += local_cell.get_total_energy();
                    topo_clusts_list.at(local_cell.get_label() - 1).add_cell(local_cell);
                }
            }
        }
    }

    int size_shares = shares.size();
    for (int ishar = 0; ishar < size_shares; ishar++)
    {
        shares.at(ishar)->set_local_max_weights(topo_clusts_list.at(shares.at(ishar)->get_1st_local_max_label() - 1).xyz_com,
                                                topo_clusts_list.at(shares.at(ishar)->get_1st_local_max_label() - 1).total_energy,
                                                topo_clusts_list.at(shares.at(ishar)->get_2nd_local_max_label() - 1).xyz_com,
                                                topo_clusts_list.at(shares.at(ishar)->get_2nd_local_max_label() - 1).total_energy);
        energy_ += shares.at(ishar)->get_total_energy();
        shares.at(ishar)->set_label(shares.at(ishar)->get_1st_local_max_label());

        topo_clusts_list.at(shares.at(ishar)->get_1st_local_max_label() - 1).add_cell(*shares.at(ishar));
        topo_clusts_list.at(shares.at(ishar)->get_2nd_local_max_label() - 1).add_cell(*shares.at(ishar));
    }

    energy_ = 0;
    for (int itopo = 0; itopo < num_of_labels; itopo++)
    {
        energy_ += topo_clusts_list.at(itopo).total_energy;
    }
}
