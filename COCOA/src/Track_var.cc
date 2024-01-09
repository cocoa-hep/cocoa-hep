#include "Track_var.hh"
#include "TRandom.h"

Track_struct::Track_struct(){};

Track_struct::Track_struct(const Track_struct &orig) = default;

Track_struct::Track_struct(int Pdgcode, int NFinalStateParticles, double Energy,
                           double Mass, double Charge,
                           double Px, double Py, double Pz,
                           double InitX, double InitY, double InitZ, bool _is_conversion_track) :
    m_lhed( -1 )
{
    is_in_endcap = false;
    is_conversion_track = _is_conversion_track;
    pdgcode = Pdgcode;
    nFinal_State_Particles = NFinalStateParticles;
    energy = Energy;
    mass = Mass;
    px = Px;
    py = Py;
    pz = Pz;
    pt = sqrt(px * px + py * py);
    absmom = sqrt(pt * pt + pz * pz);
    charge = Charge;
    initX = InitX;
    initY = InitY;
    initZ = InitZ;
    Is_track_reconstracted = true;
    Is_inside_R = true;
    Is_reach_calorimeter = true;

    eta.clear();
    phi.clear();
    x_mid_layer.clear();
    y_mid_layer.clear();
    z_mid_layer.clear();
    ind_eta.clear();
    ind_phi.clear();
    index_of_closest_topoclusters.clear();
    Rprime_to_closest_topoclusters.clear();

    rho = 0;
    IsProductInsideRadius();
}
void Track_struct::smearing()
{
    CSVReader &csv = CSVReader::GetInstance();
    // auto runData
    // 	= static_cast<DataStorage*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    float sigma_a0 = csv.hits_sigma_A0.GetBinContent(csv.hits_sigma_A0.FindBin(pt));
    float sigma_z0 = csv.hits_sigma_Z0.GetBinContent(csv.hits_sigma_Z0.FindBin(pt));
    float sigma_qp = csv.hits_sigma_QP.GetBinContent(csv.hits_sigma_QP.FindBin(pt));
    float sigma_theta = csv.hits_sigma_Theta.GetBinContent(csv.hits_sigma_Theta.FindBin(pt));
    float sigma_phi = csv.hits_sigma_Phi.GetBinContent(csv.hits_sigma_Phi.FindBin(pt));
    a0 += sigma_a0 * gRandom->Gaus(0, 1);
    z0 += sigma_z0 * gRandom->Gaus(0, 1);
    q_p += sigma_qp * gRandom->Gaus(0, 1);
    theta += sigma_theta * gRandom->Gaus(0, 1);
    phiHelix += sigma_phi * gRandom->Gaus(0, 1);
    phiHelix = phiHelix > M_PI ? phiHelix - 2 * M_PI : phiHelix < -M_PI ? phiHelix + 2 * M_PI : phiHelix;
}
void Track_struct::IsProductInsideRadius()
{
    float Rtrack = sqrt(sqr(initX) + sqr(initY));
    if (is_conversion_track){
        //Allow for conversion tracks to originate as far out as between Str0 and Str1
        if  ((Rtrack > r_inn_trkStr1) || (fabs(initZ) > pos_EndCap_trkStr1))
            Is_inside_R = false;
    }
    else if ((Rtrack > r_inn_trkPix1) || (fabs(initZ) > pos_EndCap_trkPix1))
        Is_inside_R = false;
}
void Track_struct::IsTrackReconstructed()
{
    CSVReader &csv = CSVReader::GetInstance();
    // auto runData
    // 	= static_cast<DataStorage*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    Is_track_reconstracted = true;
    float randnumb = gRandom->Uniform(1);

    float effreco = csv.hist_recon_eff.GetBinContent(csv.hist_recon_eff.FindBin(pt));
    if (randnumb < effreco)
        Is_track_reconstracted = true;
    else
        Is_track_reconstracted = false;
}

bool Track_struct::Is_Track_Useable()
{
    return Is_track_reconstracted && Is_inside_R && Is_reach_calorimeter;
}
