Event Description
====================================

The information in each event is stored as a row in a TTree object in the output ROOT file, with the branches documented in the following.


Cells
------------------------------------

**Note :** only cells contained in topoclusters are stored in the output.

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

**cell_parent_idx :** index of parent particle which contributed the largest fraction of this cell's energy

**cell_conv_el_idx :** index of conversion electron which contributed energy to this cell

**cell_parent_list :** list of parent particles which contributed to this cell's energy in decreasing order of contribution

**cell_parent_energy :** list of energies that parent particles contributed to this cell in decreasing order of contribution

Tracks
------------------------------------

The tracks are built using the Geant truth information and receive efficiency corrections and smearing of the reconstructed track parameters as specified in the configuration files in COCOA/tracking_configuration

**track_parent_idx :** index of parent particle corresponding to this track

**track_pdgid :** the PDG ID of the charged particle which creates the track

**track_d0 :** track perigee parameter d0, distance of perigee point from z-axis

**track_z0 :** track perigee parameter z0, distance of perigee point from xy-plane

**track_theta :** track perigee parameter theta

**track_phi :** track perigee parameter phi

**track_qoverp :** track perigee parameter q/p

**track_reconstructed :** track passes reconstruction efficiency cut

**track_in_acceptance :** track is inside the detector acceptance (reaches calorimeter and production vertex is inside ITS pixel layer 1)

**track_x_layer_i :** x-cordinate of the extrapolated track in the i-th calorimeter-layer

**track_y_layer_i :** y-cordinate of the extrapolated track in the i-th calorimeter-layer

**track_z_layer_i :** z-cordinate of the extrapolated track in the i-th calorimeter-layer

Nodes (of truth particle record)
------------------------------------

**node_{} :** values for constructing the graph describing the truth particle record

Truth particles
------------------------------------

**particle_pdgid :** vector storing PDG-ID of all the stable truth-particles in the event

**particle_isIso :** vector storing if the particle is isolated. If true, the value is 1. 

**particle_pt :** transverse momentum of the stable particles

**particle_eta :** pseudo-rapidity of the stable particles

**particle_phi :** azimuthal angle of the stable particles

**particle_e :** energy of the stable particles

**particle_prod_x :** x-coordinate of the production vertex of the stable particles

**particle_prod_y :** y-coordinate of the production vertex of the stable particles

**particle_prod_z :** z-coordinate of the production vertex of the stable particles

**particle_track_idx :** index of the track corresponding to this particle

**particle_eta_extrap_calo :** pseudo-rapidity value at the point where the particle enters the calorimeter

**particle_phi_extrap_calo :** azimuthal angle at the point where the particle enters the calorimeter

**particle_eta_extrap_its :** same as above, for the end of the ITS

**particle_phi_extrap_its :** same as above, for the end of the ITS

**particle_dep_energy :** amount of energy deposited by this particle in the detector

Conversion electrons
------------------------------------

**conv_el_primary_photon_idx :** index of the primary photon which conversion produced this electron

**conv_el_q :** charge of the conversion electron

**conv_el_p{} :** momentum in {x,y,z} direction of the conversion electron

**conv_el_prod_{} :** {x,y,z}-coordinate of the conversion vertex associated with this electron

Topoclusters
------------------------------------

**topo_idx :** index of topocluster (starts at 1)

**topo_bary_eta :** pseudo-rapidity of the topocluster barycenter

**topo_bary_phi :** azimuthal angle of the topocluster barycenter

**topo_bary_rho :** radius in transverse plane of the topocluster barycenter

**topo_bary_sigma_eta :** width of topocluster in eta dimension

**topo_bary_sigma_phi :** width of topocluster in phi dimension

**topo_e :** total energy of topocluster

Superclusters
------------------------------------

**supercluster_e :** total energy of supercluster

**supercluster_phi :** azimuthal angle of supercluster

**supercluster_eta :** pseudo-rapidity of supercluster

**supercluster_seed_e :** energy of topocluster used as seed for the supercluster

**supercluster_track_{} :** properties of the track associated with the seed used for the supercluster

**supercluster_topos :** list of topoclusters included in this supercluster

**supercluster_track :** index of track associated with supercluster seed

**supercluster_conv_track :** index of conversion track in case of supercluster for photon conversion

**supercluster_pdgid :** either 22 or +/-11 depending on the truth link of the track or the conversion vertex

**supercluster_N :** number of clusters included in this supercluster

Graph edges
------------------------------------

The calculation of the graph edges is defined in the files in COCOA/config/

**track_to_cell_edge_start :** list of track indices used as source nodes in track-to-cell edges

**track_to_cell_edge_end :** list of cell indices used as destination nodes in track-to-cell edges

**cell_to_cell_edge_start :** list of cell indices used as source nodes in cell-to-cell edges

**cell_to_cell_edge_end :** list of cell indices used as destination nodes in cell-to-cell edges

**particle_to_node_idx :** list of truth links to node j which particle i contributed energy to

**particle_to_node_weight :** list of the fractions of node j's energy which particle i contributed

Jets
------------------------------------

The clustering parameters of the jets are defined in the files in COCOA/config/

**true_jet_pt :** transverse component of momentum 3-vector of the jet built from truth particle

**true_jet_eta :** pseudo-rapidity of momentum 3-vector of the jet built from truth particles

**true_jet_phi :** azimuthal angle of momentum 3-vector of the jet built from truth particles

**true_jet_m :** mass of momentum 3-vector of the jet built from truth particles

**topo_jet_pt :** transverse component of momentum 3-vector of the jet built from topoclusters

**topo_jet_eta :** pseudo-rapidity of momentum 3-vector of the jet built from topoclusters

**topo_jet_phi :** azimuthal angle of momentum 3-vector of the jet built from topoclusters

**topo_jet_m :** mass of momentum 3-vector of the jet built from topoclusters

**topo_jet_constituents_list :** list of topoclusters comprising the topo jet
