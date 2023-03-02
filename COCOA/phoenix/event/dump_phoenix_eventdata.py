import json
import uproot
import numpy as np
import pickle as pkl
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse
import seaborn as sns

cmap_topos = sns.color_palette("Set3", as_cmap=True)

class TrackHelix(object):
    
    def __init__(self,d0,z0,theta,phi0,qoverp,vertex):
        
        self.theta = theta
        self.phi0 = phi0
        self.d0 = d0
        self.z0 = d0
        self.vx = vertex[0]
        self.vy = vertex[1]
        self.vz = vertex[2]
        self.qoverp = qoverp
        self.q      = np.sign(qoverp)
        self.Bz = 3.8 #0.000001
        self.rho = ( (np.sin(self.theta))/(self.qoverp*self.Bz) )*(1.0/0.299792)
        
        if d0==-10000000 or z0==-10000000:
            d0_new , z0_new = self.get_d0_z0()
            self.d0 = d0_new if d0==-10000000 else d0
            self.z0 = z0_new if z0==-10000000 else z0

    def get_d0_z0(self):
        x0 = self.vx + self.q * self.rho * (np.sin(self.phi0))
        y0 = self.vx - self.q * self.rho * (np.cos(self.phi0))
        r0 = np.sqrt(x0*x0 + y0*y0)
        a0sqr = (r0 - self.rho)**2

        # Perigee parameters
        d0 = self.q * (np.sqrt(a0sqr))
        if (self.vx**2 + self.vy**2 - a0sqr) < 0:  # numerical stability
            a0sqr = 0
        z0 = self.vz - 2 * self.rho * (1/np.tan(self.theta)) * np.arcsin(np.sqrt((self.vx**2 + self.vy**2 - a0sqr) / (4 * (self.rho**2) + 4 * self.rho * self.q * d0)))

        return d0,z0

    #parametrisation of helix, starting from perigee point
    def x_of_phi(self,phi):
        return self.vx+self.d0*np.cos(self.phi0+np.pi/2.0)+self.rho*(np.cos(phi+np.pi/2)-np.cos(self.phi0+np.pi/2))
    def y_of_phi(self,phi):
        return self.vy+self.d0*np.sin(self.phi0+np.pi/2.0)+self.rho*(np.sin(phi+np.pi/2)-np.sin(self.phi0+np.pi/2))
    def z_of_phi(self,phi):
        return self.vz+self.z0-(self.rho)*(1.0/np.tan(self.theta))*(phi-self.phi0)

def get_track_traj(d0, z0, theta, phi0, qoverp, vertex=[0,0,0]):

    track = TrackHelix(d0,z0,theta,phi0,qoverp,vertex)

    phis = np.linspace(phi0,phi0-np.sign(qoverp)*np.pi*1.0,550)
    
    xs = np.array([track.x_of_phi(phi) for phi in phis])
    ys = np.array([track.y_of_phi(phi) for phi in phis])
    zs = np.array([track.z_of_phi(phi) for phi in phis])
    
    rs = np.sqrt( xs**2+ys**2 )

    stop_idx = np.where(np.logical_or(rs >1500,np.abs(zs)>3193.9))[0]

    if len(stop_idx) > 0:
        stop_idx = stop_idx[0]

    else:
        return []

    traj = np.column_stack([xs[:stop_idx],ys[:stop_idx],zs[:stop_idx]])
    traj = traj.tolist()
    
    return traj


#def get_part_traj(px,py,pz,x0,y0,z0):
def get_part_traj(pt,eta,phi,x0,y0,z0):

    #ptot  = np.sqrt( px * px + py * py + pz * pz )
    #theta = 2*np.atan(np.e^(-1*eta)) #np.arccos( pz / ptot )
    #phi   = np.arctan2(py, px)

    ptot = pt * np.cosh(eta);
    pz   = pt * np.sinh(eta);
    px   = pt * np.cos(phi);
    py   = pt * np.sin(phi);

    times = range(10000)

    xs = np.array([x0 + t*px/ptot for t in times])
    ys = np.array([y0 + t*py/ptot for t in times])
    zs = np.array([z0 + t*pz/ptot for t in times])
    
    rs = np.sqrt( xs**2+ys**2 )

    stop_idx = np.where(np.logical_or(rs >1500,np.abs(zs)>3193.9))[0]

    if len(stop_idx) > 0:
        stop_idx = stop_idx[0]

    else:
        return []

    traj = np.column_stack([xs[:stop_idx],ys[:stop_idx],zs[:stop_idx]])
    traj = traj.tolist()
    
    return traj

def get_topoclusters(cell_topo_idx, cell_eta, cell_phi, cell_e):
    
    if len(cell_topo_idx) == 0:
        return []
    n_topoclusters = max(cell_topo_idx)

    topoclusters = []
    for i in range(n_topoclusters):
        mask = cell_topo_idx == i+1
        
        tmp_e = cell_e[mask].sum()
        tmp_eta = (cell_eta * cell_e)[mask].sum() / tmp_e
        tmp_phi = (cell_phi * cell_e)[mask].sum() / tmp_e

        tmp_tc = {
            'energy': tmp_e.item(),
            'eta'   : tmp_eta.item(),
            'phi'   : tmp_phi.item()
        }
        
        topoclusters.append(tmp_tc)
        
    return topoclusters

def get_cells(cell_x, cell_y, cell_z, cell_e):
    
    cells = []
    
    for i,x in enumerate(cell_x):
        size = np.log(cell_e[i].item())*5
        cells.append({
            'type': 'Point',
            'pos': [cell_x[i].item(),cell_y[i].item(),cell_z[i].item()], #,size,size,size],
            #'x': cell_x[i].item(),
            #'y': cell_y[i].item(),
            #'z': cell_z[i].item(),
            'color': '#FFFF00', #yellow
        })
        
    return cells

def energy_transparency(cell_e):
    
    return (np.log(cell_e) - np.min(np.log(cell_e)))/np.ptp(np.log(cell_e))

def get_color(layer,scaleby,topoidx=-1):

    #lay_max = np.array([
    #                    [0, 255, 128],
    #                    [0, 255, 128],
    #                    [0, 255, 128],
    #                    [127, 0, 255],
    #                    [127, 0, 255],
    #                    [127, 0, 255]
    #                   ])
    
    lay_max = np.array([[110, 222, 138],
                        [146, 230, 167],
                        [183, 239, 197],
                        [129, 137, 255],
                        [150, 162, 255],
                        [174, 184, 255]
                       ])

    if(topoidx > 0):
        topoidx = topoidx + 1 # avoid first (yellowey) color in cmap
        while topoidx>12:
            topoidx = topoidx - 12
        rgb = [int(cmap_topos(topoidx)[i]*255) for i in range(3)]
    else:
        rgb = lay_max[layer]
        rgb = rgb*scaleby
        rgb = np.clip(rgb,0,255).astype(int).tolist()

    #return f"rgb({rgb[0]}, {rgb[1]}, {rgb[2]})"
    return rgb

def get_vertices(cell_x, cell_y, cell_z, cell_l, cell_e, cell_topoidx, vertices_df, dump_vertices=False, dump_centres=False):
    
    if len(cell_x) == 0:
        return []
    merge_supercells = True
    size = 0.2
    vertices = []
    hashlist = []
    distlist = []
    colors   = ["rgb(110, 222, 138)","rgb(146, 230, 167)","rgb(183, 239, 197)","rgb(129, 137, 255)","rgb(150, 162, 255)","rgb(174, 184, 255)"]
    layer_noise = [13, 34.,41.,75.,50.,25.]
    cell_noise = cell_l.copy()
    for layer, noise in enumerate(layer_noise):
        cell_noise[cell_l==layer] = noise

    SNR = cell_e / cell_noise
    SNR_01 = np.clip(SNR/4.6,0.25,1)
    #SNR_01 = np.clip(2*SNR/4.6,0.4,1)
    #SNR_01 = (SNR - np.min(SNR)) / np.ptp(SNR)
    #SNR_01 = np.clip(SNR_01,0.0,1.0)
    #SNR_01 = np.power(np.arctan(SNR_01*np.pi/2),1./2)
    #SNR_01 = np.clip(SNR_01,0.1,1.0)
    energy_01 = (np.log(cell_e) - np.min(np.log(cell_e)))/np.ptp(np.log(cell_e))

    event_cells = np.column_stack([cell_x,cell_y,cell_z])
    
    for idx, cell in enumerate(tqdm(event_cells)):
        layer_df = vertices_df[vertices_df['layer']==cell_l[idx]]
        if len(layer_df) == 0: continue
        mpos = np.column_stack([layer_df['vmx'].to_numpy(),layer_df['vmy'].to_numpy(),layer_df['vmz'].to_numpy()])
        pos  = np.repeat(np.expand_dims(cell,axis=0),len(mpos),axis=0)
        dist = np.linalg.norm(pos - mpos,axis=1)

        if(len(dist)==0):
            print(layer_df)
            print('cell_l=',cell_l[idx])
        i_min = np.argmin(dist)
        #if dist[i_min]>5:
        #    continue
            
        distlist.append(dist[i_min])
        winner = layer_df.iloc[i_min]['hash'].item()
        if winner in hashlist:
            print('WARNING: found duplicate! hash=',winner)
            continue
        hashlist.append(winner)
        rows   = layer_df[layer_df['hash']==winner]
        rows   = rows.sort_values('eta_idx')

    
        if merge_supercells and not dump_vertices and not dump_centres:

            flat_coords = []
            
            if len(rows)>1:
                vlists = [
                    [0,2,4,6],
                    [1,3,5,7],
                ]

                first_and_last = [0,len(rows)-1] if len(rows)>1 else [0]

                for subidx,subcell in enumerate(first_and_last):
                    row = rows.iloc[subcell]
                    for v in vlists[subidx] :
                        flat_coords = flat_coords + [row['v{}x'.format(v)].item(),row['v{}y'.format(v)].item(),row['v{}z'.format(v)].item()]
            else:
                row = rows
                for v in range(8):
                    flat_coords = flat_coords + [row['v{}x'.format(v)].item(),row['v{}y'.format(v)].item(),row['v{}z'.format(v)].item()]                
            

            layer = int(row['layer'].item())
            energy = cell_e[idx].item()
            #if layer > 2:
            #    layer = layer - 1
            vertices.append({
                'type': 'IrregularCaloCells',
                'layer': layer,
                'energy': energy,
                'vtx': flat_coords,
                'cluster': cell_topoidx[idx].item(),
                'color': get_color(layer, 1.0, cell_topoidx[idx].item()),
                'opacity': SNR_01[idx] #energy_01[idx].item(),
            })
        else:

            for _, row in rows.iterrows():

                if dump_centres: #dump the geometric centre of each cell as a point to display as a hit
                    v = 'm'
                    vertices.append({
                        'type': 'Point',
                        'pos': [row['v{}x'.format(v)].item(),row['v{}y'.format(v)].item(),row['v{}z'.format(v)].item()],
                        'color': '#FFFF00', #yellow
                    })
                elif dump_vertices: #dump each vertex as a separate point to display as a hit
                    for v in range(8):
                        vertices.append({
                            'type': 'Point',
                            'pos': [row['v{}x'.format(v)].item(),row['v{}y'.format(v)].item(),row['v{}z'.format(v)].item()],
                            'color': '#FFFF00', #yellow
                        })
                else: #dump full coords per cell to display as a polyhedron
                    flat_coords = []
                    vlist = range(8)
                    for v in vlist:
                        flat_coords = flat_coords + [row['v{}x'.format(v)].item(),row['v{}y'.format(v)].item(),row['v{}z'.format(v)].item()]

                    layer = int(row['layer'].item())
                    energy = cell_e[idx].item()
                    #if layer > 2:
                    #    layer = layer - 1
                    vertices.append({
                        'type': 'IrregularCaloCells',
                        'layer': layer,
                        'energy': energy,
                        'vtx': flat_coords,
                        'color': colors[layer],
                        'opacity': 1.0 #energy_01[idx].item(),
                    })

    return vertices

def get_scaled_jet_e(e):
    x_, y_ = 20_000, 300_000
    x,  y  = 25_000,  40_000
    
    scale = (y - x) / (y_ - x_)
    scaled_e = (e - x_) * scale + x
    
    return scaled_e

def cocoa_to_phoenix(ntuple_path, cell_path, output_path, nevents=-1, firstevent=0):
    
    tree = uproot.open(ntuple_path)['Out_Tree']

    if nevents<0:
        nevents = tree.num_entries
    lastevent = min(tree.num_entries, firstevent + nevents)

    # truth particles
    part_pdgid = tree["particle_pdgid"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_pt = tree["particle_pt"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_eta = tree["particle_eta"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_phi = tree["particle_phi"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    # part_px = tree["particle_px"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    # part_py = tree["particle_py"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    # part_pz = tree["particle_pz"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_prod_x = tree["particle_prod_x"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_prod_y = tree["particle_prod_y"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    part_prod_z = tree["particle_prod_z"].array(library='np', entry_stop=lastevent,entry_start=firstevent)

    # tracks
    track_d0 = tree["track_d0"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    track_z0 = tree["track_z0"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    track_theta = tree["track_theta"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    track_phi = tree["track_phi"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    track_qoverp = tree["track_qoverp"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    track_pdgid = tree["track_pdgid"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    

    # photon conversion tracks
    conv_el_q = tree["conv_el_q"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_px = tree["conv_el_px"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_py = tree["conv_el_py"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_pz = tree["conv_el_pz"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_prod_x = tree["conv_el_prod_x"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_prod_y = tree["conv_el_prod_y"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    conv_el_prod_z = tree["conv_el_prod_z"].array(library='np', entry_stop=lastevent,entry_start=firstevent)

    # topoclusters
    cell_topo_idx = tree["cell_topo_idx"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    
    #cells
    cell_x   = tree["cell_x"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    cell_y   = tree["cell_y"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    cell_z   = tree["cell_z"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    cell_l   = tree["cell_layer"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    cell_eta = tree["cell_eta"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    cell_phi = tree["cell_phi"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    cell_e = tree["cell_e"].array(library='np', entry_stop=lastevent,entry_start=firstevent)    
    
    # jets
    true_jet_pt  = tree["true_jet_pt"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    true_jet_eta = tree["true_jet_eta"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    true_jet_phi = tree["true_jet_phi"].array(library='np', entry_stop=lastevent,entry_start=firstevent)
    true_jet_m   = tree["true_jet_m"].array(library='np', entry_stop=lastevent,entry_start=firstevent)

    ### Retrieve map of cells
    cell_df   = pkl.load(open(cell_path,'rb'))

    ### Make up for iron gap <-- already converted upstream
    cell_df.loc[ cell_df['layer']==4, 'layer'] = 3
    cell_df.loc[ cell_df['layer']==5, 'layer'] = 4
    cell_df.loc[ cell_df['layer']==6, 'layer'] = 5

    ### Calculate geometric means of each cell from its 8 vertices
    cell_df['vmx'] = np.mean(np.column_stack([cell_df['v{}x'.format(i)].to_numpy() for i in range(8)]),axis=1)
    cell_df['vmy'] = np.mean(np.column_stack([cell_df['v{}y'.format(i)].to_numpy() for i in range(8)]),axis=1)
    cell_df['vmz'] = np.mean(np.column_stack([cell_df['v{}z'.format(i)].to_numpy() for i in range(8)]),axis=1)

    event_data = {}
    
    for i in range(nevents):
                
        event_dict = dict()
        
        event_dict['event number']: i
        event_dict['run number']  : 0

            
        #    
        # tracks and trajectories
        #------------------
        
        n_tracks = track_d0[i].shape[0]
        
        tracks = {}
        all_tracks = []
        el_tracks = []; chhad_tracks = []; mu_tracks = []

        for j in range(n_tracks):
            track_traj = get_track_traj(
                    track_d0[i][j], track_z0[i][j], track_theta[i][j], track_phi[i][j], track_qoverp[i][j])

            tmp_track = {
                'pos': track_traj,
                'color': '0xff0000'
            }
            all_tracks.append(tmp_track)
            
            if track_pdgid[i][j] in [11, -11] :
                tmp_track = {
                    'pos': track_traj,
                    'color': '0xff8700'
                }
                el_tracks.append(tmp_track)

            elif track_pdgid[i][j] in [13, -13] :
                tmp_track = {
                    'pos': track_traj,
                    'color': '0xffd100'
                }
                mu_tracks.append(tmp_track)

            else:
                tmp_track = {
                    'pos': track_traj,
                    'color': '0xff5d00'
                }
                chhad_tracks.append(tmp_track)

        # Add displaced tracks of electrons from photon conversions
        n_conv_el = conv_el_prod_x[i].shape[0]
        
        conv_el_tracks = []

        for j in range(n_conv_el):

            #conv_el_d0, conv_el_z0 = get_d0_z0_from_xyz(conv_el_prod_x[i][j] ,conv_el_prod_y[i][j],conv_el_prod_z[i][j],conv_el_px[i][j],conv_el_py[i][j],conv_el_pz[i][j],1,3.8) #TODO charge
            conv_el_p = np.sqrt(conv_el_px[i][j]**2 + conv_el_py[i][j]**2 + conv_el_pz[i][j]**2)
            conv_el_theta = np.arccos(conv_el_pz[i][j]/conv_el_p)
            conv_el_phi   = np.arctan2(conv_el_py[i][j],conv_el_px[i][j])
            conv_el_traj = get_track_traj(
                    0, 0, conv_el_theta, conv_el_phi, conv_el_q[i][j]/conv_el_p, [conv_el_prod_x[i][j],conv_el_prod_y[i][j],conv_el_prod_z[i][j]])

            tmp_track = {
                'pos': conv_el_traj,
                'color': '0xff0000'
            }
            conv_el_tracks.append(tmp_track)

        # Add neutral particle trajectories as "tracks" with no curvature and drawn dashed (TODO)
        ph_trajs = []; nuhad_trajs = [];

        for j in range(part_pdgid[i].shape[0]):
            if part_pdgid[i][j] == 22 or part_pdgid[i][j] == 111 :

                #part_traj = get_part_traj(part_px[i][j], part_py[i][j], part_pz[i][j], part_prod_x[i][j], part_prod_y[i][j], part_prod_z[i][j])
                part_traj = get_part_traj(part_pt[i][j], part_eta[i][j], part_phi[i][j], part_prod_x[i][j], part_prod_y[i][j], part_prod_z[i][j])

                tmp_traj = {
                    'pos': part_traj,
                    'color': '0xff0000'
                }
                ph_trajs.append(tmp_traj)

        tracks['all_tracks']   = all_tracks
        tracks['el_tracks']     = el_tracks
        tracks['chhad_tracks'] = chhad_tracks
        tracks['mu_tracks']    = mu_tracks
        tracks['conv_el_tracks']= conv_el_tracks
        tracks['ph_trajs']     = ph_trajs

        event_dict['Tracks'] = tracks
        
                
        #    
        # jets
        #------------------           

        jets = {}

        n_true_jets = true_jet_pt[i].shape[0]
        true_jets = []

        for j in range(n_true_jets):
            jet_p = true_jet_pt[i][j] * np.cosh(true_jet_eta[i][j])   
            jet_e = np.sqrt(true_jet_m[i][j] **2 + jet_p **2)
            
            tmp_jet = {
                'eta': true_jet_eta[i][j].item(),
                'phi': true_jet_phi[i][j].item(),
                'coneR': 0.4,
                'energy': get_scaled_jet_e(jet_e.item()),
#                     'color': '0x72fcfc'
            }
            true_jets.append(tmp_jet)
        jets['truth jets'] = true_jets
            
        
        event_dict['Jets'] = jets
        
        
                
        #    
        # calclusters
        #------------------
        
        calo_clusters = {}
        calo_clusters['topoclusters'] = get_topoclusters(cell_topo_idx[i], cell_eta[i], cell_phi[i], cell_e[i])        
        event_dict['CaloClusters'] = calo_clusters

        #
        # hits
        #------------------
        
        cells = {}
        cells['hits'] = get_cells(cell_x[i],cell_y[i],cell_z[i],cell_e[i])
        event_dict['Hits'] = cells
        
        #
        # cells
        #------------------
        
        cells = {}
        cells['centres'] = get_vertices(cell_x[i],cell_y[i],cell_z[i],cell_l[i],cell_e[i],cell_topo_idx[i],cell_df,dump_centres=True)
        event_dict['Hits'] = cells

        
        #
        # COCOACaloCells
        #------------------
        
        scc = {}
        scc['vertices'] = get_vertices(cell_x[i],cell_y[i],cell_z[i],cell_l[i],cell_e[i],cell_topo_idx[i],cell_df)
        event_dict['IrregularCaloCells'] = scc
        
        
        # event
        event_data[f'event_num_{i}'] = event_dict
        
    json_object = json.dumps(event_data, indent=2)
    with open(output_path, "w") as outfile:
        outfile.write(json_object)

    print("Results written to ",output_path)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",   dest="input",   type=str, help="path to input ROOT file", required=True)
    parser.add_argument("-o","--output",  dest="output",  type=str, help="path to output json file", default="events.json")
    parser.add_argument("-c","--cells",   dest="cells",   type=str, help="path to cell geometry lookup table (pkl)", default="cells_cocoa_default.pkl")
    parser.add_argument("-n","--nevents", dest="nevents", type=int, help="number of events to parse", default=1)
    args = parser.parse_args()

    cocoa_to_phoenix(args.input, args.cells, args.output, args.nevents)


if __name__ == '__main__':
    main()
