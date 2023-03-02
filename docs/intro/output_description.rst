Event Description
====================================

The event information is stored in the output ROOT file. The cells within a topocluster are only stored.

The output ROOT file has the following branches, containing informations at per event level.
------------------------------------

**cell_layer :** Which calorimeter layer the cell belongs to. 0 corresponds to ECAL1 and 5 corresponds to HCAL3

**cell_x :** x-coordinate of the cell

**cell_y :** y-coordinate of the cell

**cell_z :** z-coordinate of the cell

**cell_eta :** pseudo-rapidity of the cell

**cell_phi :** azimuthal angle of the cell

**cell_e :** total energy deposited in the cell

**cell_che :** total energy deposited in the cell by calorimeter shower component, which was initiated by a charged particle

**cell_nue :** total energy deposited in the cell by calorimeter shower component, which was initiated by a neutral particle

**cell_topo_idx :** Topocluster label to which the cell belongs

**track_pflow_object_idx :** the index of pflow object to a given track

**track_pdgid :** the PDG ID of the charged particle which creates the track

**track_d0 :** track perigee parameter d0, distance of perigee point from z-axis

**track_z0 :** track perigee parameter z0, distance of perigee point from xy-plane

**track_theta :** track perigee parameter theta

**track_phi :** track perigee parameter phi

**track_qoverp :** track perigee parameter q/p

**track_eta_layer_i :** pseudo-rapidity of the extrapolated track in the i-th calorimeter-layer

**track_phi_layer_i :** azimuthal angle of the extrapolated track in the i-th calorimeter-layer

**track_x_layer_i :** x-cordinate of the extrapolated track in the i-th calorimeter-layer

**track_y_layer_i :** y-cordinate of the extrapolated track in the i-th calorimeter-layer

**track_z_layer_i :** z-cordinate of the extrapolated track in the i-th calorimeter-layer

**particle_pdgid :** vector storing PDG-ID of all the stable truth-particles in the event

**particle_isIso :** vector storing if the particle is isolated. If true, the value is 1. 

**particle_px :** x-component of the momentum 3-vector of the stable particles

**particle_py :** y-component of the momentum 3-vector of the stable particles

**particle_pz :** z-component of the momentum 3-vector of the stable particles

**particle_e :** energy of the stable particles

**particle_prod_x :** x-coordinate of the production vertex of the stable particles

**particle_prod_y :** y-coordinate of the production vertex of the stable particles

**particle_prod_z :** z-coordinate of the production vertex of the stable particles

**particle_eta_layi :** extrapolated pseudo-rapidity of the particle in the i-th calorimeter layer

**particle_phi_layi :** extrapolated azimuthal angle of the particle in the i-th calorimeter layer

**pflow_eta :** pseudo-rapidity of the particle-flow object

**pflow_phi :** azimuthal angle of the particle-flow object

**pflow_px :** x-component of the momentum 3-vector of the particle-flow object

**pflow_py :** y-component of the momentum 3-vector of the particle-flow object

**pflow_pz :** z-component of the momentum 3-vector of the particle-flow object

**pflow_e :** energy of the particle-flow object

**pflow_charge :** electric charge of the particle-flow object

**pflow_jet_pt :** transverse component of momentum 3-vector of the jet built from pflow constituents

**pflow_jet_eta :** pseudo-rapidity of momentum 3-vector of the jet built from pflow constituents

**pflow_jet_phi :** azimuthal angle of momentum 3-vector of the jet built from pflow constituents

**pflow_jet_m :** mass of momentum 3-vector of the jet built from pflow constituents

**topo_jet_pt :** transverse component of momentum 3-vector of the jet built from topoclusters

**topo_jet_eta :** pseudo-rapidity of momentum 3-vector of the jet built from topoclusters

**topo_jet_phi :** azimuthal angle of momentum 3-vector of the jet built from topoclusters

**topo_jet_m :** mass of momentum 3-vector of the jet built from topoclusters




