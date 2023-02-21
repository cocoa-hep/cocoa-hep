#include "Supercluster.hh"
#include "Topo_clust_var.hh"
#include "Superclustering.hh"


ConversionVertex::ConversionVertex()
{
    track_list = {};
    topo_list = {};
    x = 0;
    y = 0;
    z = 0;
}

bool ConversionVertex::contains(Topo_clust topo)
{
    return (std::find(topo_list.begin(),topo_list.end(),topo.label) != topo_list.end());
}

bool ConversionVertex::contains(Track_struct track)
{
    return (std::find(track_list.begin(),track_list.end(),track.position_in_list) != track_list.end());
}

bool ConversionVertex::try_add_track(Track_struct track)
{
    bool success = false;
    if(track.is_conversion_track)
    {
        if(track_list.size()==0)
        {
            track_list.push_back(track.position_in_list);
            x = track.initX;
            y = track.initY;
            z = track.initZ;
            success = true; //first track added
        }
        else if(contains(track))
        {
            success = false; //already is there
        }
        else if(track_list.size()==1)
        {
            float dist2 = sqr(track.initX - x)
                        + sqr(track.initY - y)
                        + sqr(track.initZ - z);
            if(sqrtf(dist2)<0.0001)
            {
                track_list.push_back(track.position_in_list);
                success = true; //second track added
            }
        }
    }

    return success;
}

bool ConversionVertex::try_add_topo(Topo_clust topo)
{
    if (topo.closest_tracks.size()==0)
        return false;

    int success = false;
    int nearest_track_idx = topo.closest_tracks.at(0).second;

    if(track_list.size() < 2)
    {
        success = false;
    }
    else if(contains(topo))
    {
        success = false;
    }
    else if(std::find(track_list.begin(),track_list.end(),nearest_track_idx) != track_list.end())
    {
        topo_list.push_back(topo.label);
        success = true;
    }
    return success;
}

int ConversionVertex::get_friend(Track_struct track)
{
    if(track_list.size() < 2 || !track.is_conversion_track)
        return -1;

    auto iter = std::find(track_list.begin(),track_list.end(),track.position_in_list);
    if(iter == track_list.end())
        return -1;

    int where = std::distance(track_list.begin(), iter);
    return track_list.at(1 - where); // return the other one (at(1) if 0, at(0) if 1)
}

int ConversionVertex::get_friend(Topo_clust topo)
{
    if(topo_list.size() < 2)
        return -1;

    auto iter = std::find(topo_list.begin(),topo_list.end(),topo.label);
    if(iter == topo_list.end())
        return -1;

    int where = std::distance(topo_list.begin(), iter);
    return topo_list.at(1 - where); // return the other one (at(1) if 0, at(0) if 1)
}


Supercluster::Supercluster()
{
    com = {0,0,0};
    eta = -999;
    phi = -999;
}

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
        total_energy += topo_members.at(itopo).EM_energy;
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

void Supercluster::set_conv_track(Track_struct track)
{
    track_conv = track;
}

void Supercluster::set_conv_vertex(ConversionVertex vtx)
{
    conv_vertex = vtx;
}

std::vector<Topo_clust> Supercluster::get_clusters()
{
    return topo_members;
}

Track_struct Supercluster::get_track()
{
    return track_assoc;
}

Track_struct Supercluster::get_conv_track()
{
    return track_conv;
}

ConversionVertex Supercluster::get_conv_vertex()
{
    return conv_vertex;
}