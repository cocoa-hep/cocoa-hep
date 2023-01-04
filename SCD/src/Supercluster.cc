#include "Supercluster.hh"
#include "Topo_clust_var.hh"


Supercluster::Supercluster(){}

void Supercluster::sort_by_energy()
{
    sort(topo_members.begin(),topo_members.end(),compare_e);
}

void Supercluster::update()
{
    n_clusters = topo_members.size();
    total_energy = 0;
    for (int itopo = 0; itopo < n_clusters; itopo++)
    {
        total_energy += topo_members.at(itopo).total_energy;
    }
}

void Supercluster::add_cluster(Topo_clust topo)
{
    topo_members.push_back(topo);
    update();
}

void Supercluster::set_seed(Topo_clust seed)
{
    if(topo_members.size() > 0) {
        topo_members.at(0) = seed;
        update();
    }
    else add_cluster(seed);
}

void Supercluster::set_track(Track_struct track)
{
    track_assoc = track;
}

Topo_clust Supercluster::get_seed()
{
    return topo_members[0];
}