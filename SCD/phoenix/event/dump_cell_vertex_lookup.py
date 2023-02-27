import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm, rc
import pandas as pd
from copy import copy, deepcopy
from tqdm import tqdm

font = {'family' : 'DejaVu Sans',
#        'weight' : 'bold',
        'size'   : 8 }

rc('font', **font)
rc('axes', linewidth=0.5)
rc('lines', lw=0.5)
rc('axes', axisbelow=False)
plt.rcParams['xtick.major.size'] = 1.0
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.size'] = 1.0
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.major.size'] = 1.0
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.size'] = 1.0
plt.rcParams['ytick.minor.width'] = 0.5

class Cell:

    def __init__(self):
        self.granularity = [256,256,128,-1,64,64,32]
        self.vertices = [[],[],[],[],[],[],[],[]]
        self.vertices_dict = {}
        self.region = 0 #Barrel (2) or endcap (-1,1)
        self.layer  = -1
        self.eta_idx = 0
        self.phi_idx = 0
        self.dphi = 0

    def RotatePhi(self):

        self.dphi = 2*np.pi/self.granularity[self.layer]
        rotated_cell = deepcopy(self)
        rotated_cell.phi_idx = self.phi_idx + 1 #increment phi division

        rotated_cell.vertices[0] = rotated_cell.vertices[4][:]
        rotated_cell.vertices[1] = rotated_cell.vertices[5][:]
        rotated_cell.vertices[2] = rotated_cell.vertices[6][:]
        rotated_cell.vertices[3] = rotated_cell.vertices[7][:]

        for i,v in enumerate(self.vertices):
            if i<4: 
                continue #only shift last 4
            x  = v[0]
            y  = v[1]
            xp = x*np.cos(self.dphi) - y*np.sin(self.dphi)
            yp = x*np.sin(self.dphi) + y*np.cos(self.dphi)
            rotated_cell.vertices[i][0] = xp
            rotated_cell.vertices[i][1] = yp

        rotated_cell = Update(rotated_cell)

        return rotated_cell

    def MirrorZ(self):

        mirrored_cell =  deepcopy(self)

        mirrored_cell.region  = int(-1*mirrored_cell.region)

        for i,v in enumerate(self.vertices):
            mirrored_cell.vertices[i] = [self.vertices[i][0], self.vertices[i][1], -1*self.vertices[i][2]]

        mirrored_cell = Update(mirrored_cell)

        return mirrored_cell


def Update(cell):

    cell.vertices_dict = {}
    for i,v in enumerate(cell.vertices):
        cell.vertices_dict['v{}x'.format(i)] = [v[0]]
        cell.vertices_dict['v{}y'.format(i)] = [v[1]]
        cell.vertices_dict['v{}z'.format(i)] = [v[2]]
    cell.vertices_dict['layer'] = cell.layer
    cell.vertices_dict['region'] = cell.region
    cell.vertices_dict['eta_idx'] = cell.eta_idx
    cell.vertices_dict['phi_idx'] = cell.phi_idx
    cell.vertices_dict['hash']    = int(np.sign(cell.region)*int( f'{abs(cell.region)}{cell.layer}{cell.eta_idx}{cell.phi_idx:03}' ))

    return cell

def EtaToTheta( eta ):
    return 2.0 * np.arctan( np.exp( -eta ) )

#--------------------------------------------------

#
# Assuming that the barrel covers half of the detector eta range
#

nBins        = 128
etaMax       = 3.0
etaMaxBarrel = etaMax / 2.0
dEta         = etaMax / nBins

linewidth = 0.5

layer_names = [ 'ECAL1', 'ECAL2', 'ECAL3', 'gap', 'HCAL1', 'HCAL2', 'HCAL3' ]

r_all       = [ 1500, 1594.11, 1970.56, 2017.62, 2097.62, 2496.17, 3585.57, 4063.84 ] # Fe gap missing

layer_mergeFactors = [ 1, 1, 2, 1, 4, 4, 8 ]

nLayers     = len( r_all ) - 1

tableau_colors = list( mcolors.TABLEAU_COLORS )
colors = tableau_colors[:3]
colors.append( "white" )
colors += tableau_colors[3:]


nBins_barrel = nBins // 2
nBins_endcap = nBins_barrel

x_min_border      = np.zeros( nBins_barrel )
x_max_border      = np.zeros( nBins_barrel )
y_low_max_border  = np.zeros( nBins_barrel )
y_high_max_border = np.zeros( nBins_barrel )

cell_idx = 0
cell_df = pd.DataFrame()

#
# Endcap
#

for iLayer in tqdm(range(nLayers)):

    eta_idx = 0
    if iLayer == 3: # Fe gap
        continue
    
    for iCell in tqdm(range(nBins_endcap)):
    
        if iCell % layer_mergeFactors[iLayer] == 0 and iCell>0:
            eta_idx = eta_idx + 1

        color = colors[iLayer]
        d     = r_all[iLayer + 1] - r_all[iLayer]
        lz    = r_all[0] / np.tan( EtaToTheta( etaMaxBarrel ) )
        
        theta_min = EtaToTheta( etaMax - iCell * dEta )
        theta_max = EtaToTheta( etaMax - ( iCell  + 1 ) * dEta )
        
        lz += ( r_all[iLayer] - r_all[0] ) * np.cos( theta_min )
        
        y_low_min = lz * np.tan( theta_min )
        y_low_max = lz * np.tan( theta_max )
    
        dz         = d * np.cos( theta_min )
        y_high_min = ( lz + dz ) * np.tan( theta_min )
        y_high_max = ( lz + dz ) * np.tan( theta_max )

        x_min_border[iCell]      = lz
        x_max_border[iCell]      = lz + dz
        y_low_max_border[iCell]  = y_low_max
        y_high_max_border[iCell] = y_high_max

        #Fill cell vertices dataframe
        cell = Cell()
        cell.layer = iLayer
        cell.region = 1
        cell.eta_idx = eta_idx
        cell.phi_idx = 0
        cell.vertices[0] = [0,y_low_min,lz]
        cell.vertices[1] = [0,y_low_max,lz]
        cell.vertices[2] = [0,y_high_min,lz+dz]
        cell.vertices[3] = [0,y_high_max,lz+dz]
        cell.vertices[4] = cell.vertices[0]
        cell.vertices[5] = cell.vertices[1]
        cell.vertices[6] = cell.vertices[2]
        cell.vertices[7] = cell.vertices[3]
        cell.dphi = -1
        cell = Update(cell)

        for slice in range(cell.granularity[cell.layer]):
            #if slice > 0: continue HACK!
            rotated_cell = cell.RotatePhi()
            rotated_cell_df_i = pd.DataFrame(rotated_cell.vertices_dict)
            if(len(cell_df)==0):
                cell_df = rotated_cell_df_i
            else:
                cell_df = pd.concat([cell_df,rotated_cell_df_i])
            cell_df = pd.concat([cell_df,pd.DataFrame(rotated_cell.MirrorZ().vertices_dict)])
            cell = rotated_cell

        if iCell == 0:
            plt.plot( [lz, lz], [y_low_min, y_low_max] , c = color, lw = linewidth, label = layer_names[iLayer] )
        else:
            plt.plot( [lz, lz], [y_low_min, y_low_max] , c = color, lw = linewidth )
        plt.plot( [lz + dz, lz + dz], [y_high_min, y_high_max], c = color, lw = linewidth )

        if iCell > 1:
            plt.plot( [lz, x_min_border[iCell - 1]], [y_low_min, y_low_max_border[iCell - 1]], c = color, lw = linewidth )
            plt.plot( [lz + dz, x_max_border[iCell - 1]], [y_high_min, y_high_max_border[iCell - 1]], c = color, lw = linewidth )
        
        if iCell == nBins_endcap - 1:
            plt.plot( [lz + dz, lz], [y_high_max, y_low_max], c = color, lw = linewidth )
        if iCell % layer_mergeFactors[iLayer] == 0:
            plt.plot( [lz, lz + dz], [y_low_min, y_high_min], c = color, lw = linewidth )
        
#
# Barrel
#

x_low_max_border  = np.zeros( nBins_barrel )
x_high_min_border = np.zeros( nBins_barrel )
x_high_max_border = np.zeros( nBins_barrel )
y_high_border     = np.zeros( nBins_barrel )
y_low_border      = np.zeros( nBins_barrel )

for iLayer in tqdm(range(nLayers)):
    
    eta_idx = 0

    for iCell in tqdm(range(nBins_barrel)):

        if iCell % layer_mergeFactors[iLayer] == 0 and iCell>0:
            eta_idx = eta_idx + 1

        color = colors[iLayer]
        r     = r_all[iLayer]
        d     = r_all[iLayer + 1] - r_all[iLayer]
    
        eta_min = iCell * dEta
        eta_max = eta_min + dEta
        theta_min = EtaToTheta( eta_min )
        theta_max = EtaToTheta( eta_max )

        if iLayer == 0:
            
            x_min_low = r / np.tan( theta_min )
            x_max_low = r / np.tan( theta_max )
            y_low     = r

        else:

            if iCell == 0:
                x_min_low = 0.0
            else:
                x_min_low = x_high_min_border[ iCell ]
            x_max_low = x_high_max_border[ iCell ]
            
            y_low = y_high_border[ iCell ]        
        
        dy     = d * np.sin( theta_min )
        dx_min = d * np.cos( theta_min )
        dx_max = dy / np.tan( theta_max )

#        if iLayer == 0:
#            print( x_max_low + dx_max )
#        else:
#            print( x_max_low )
#
#        print()

        x_low_max_border[iCell]  = x_max_low
        x_high_min_border[iCell] = x_min_low + dx_min
        x_high_max_border[iCell] = x_max_low + dx_max
        y_low_border[iCell]      = y_low
        y_high_border[iCell]     = y_low + dy
        
        if iLayer == 3: # Fe gap
            continue

        #Fill cell vertices dataframe
        cell = Cell()
        cell.layer = iLayer
        cell.region = 2
        cell.eta_idx = eta_idx
        cell.phi_idx =0 
        cell.vertices[0] = [0,y_low,x_min_low]
        cell.vertices[1] = [0,y_low,x_max_low]
        #cell.vertices[2] = [0,y_high_border[iCell - 1],x_high_max_border[iCell - 1]]
        cell.vertices[2] = [0,y_high_border[iCell],x_high_min_border[iCell]]
        cell.vertices[3] = [0,y_high_border[iCell],x_high_max_border[iCell]]
        cell.vertices[4] = cell.vertices[0]
        cell.vertices[5] = cell.vertices[1]
        cell.vertices[6] = cell.vertices[2]
        cell.vertices[7] = cell.vertices[3]
        cell.dphi = -1
        cell = Update(cell)

        for slice in range(cell.granularity[cell.layer]):
            #if slice > 0: continue HACK!
            rotated_cell = cell.RotatePhi()
            rotated_cell_df_i = pd.DataFrame(rotated_cell.vertices_dict)
            if(len(cell_df)==0):
                cell_df = rotated_cell_df_i
            else:
                cell_df = pd.concat([cell_df,rotated_cell_df_i])
            cell_df = pd.concat([cell_df,pd.DataFrame(rotated_cell.MirrorZ().vertices_dict)])
            cell = rotated_cell



        plt.plot( [ x_min_low, x_max_low ], [ y_low, y_low ], c = color, lw = linewidth )
        if ( iCell % layer_mergeFactors[iLayer] == 0 ) and iCell > 0:
            plt.plot( [ x_min_low, x_high_max_border[ iCell - 1 ] ], [ y_low, y_high_border[iCell - 1] ], c = color, lw = linewidth )
        elif iCell > 0:
#            plt.plot( [ x_min_low, x_min_low + dx_min ], [ y_low, y_low + dy ], c = "k", lw = linewidth )
            plt.plot( [ x_high_max_border[iCell - 1], x_high_min_border[iCell] ], [ y_high_border[iCell - 1], y_high_border[iCell] ], c = color, lw = linewidth )
            plt.plot( [ x_min_low, x_low_max_border[iCell - 1] ], [ y_low, y_low_border[iCell - 1] ], c = color, lw = linewidth )
        if iCell == nBins_barrel - 1:
            plt.plot( [x_max_low, x_max_low + dx_max], [y_low, y_low + dy], c = color, lw = linewidth )
            
#        plt.plot( [ x_max_low, x_max_low + dx_max ], [ y_low, y_low + dy ], c = color, lw = linewidth )
        plt.plot( [ x_min_low + dx_min, x_max_low + dx_max ], [ y_low + dy, y_low + dy ], c = color, lw = linewidth )

#
# ID
#

### dump cell vertices on plot
for v in [2,3]:
    plot_df = cell_df[ (cell_df['phi_idx']==1) & (cell_df['region']>0) ]
    plt.scatter(plot_df[f'v{v}z'],plot_df[f'v{v}y'])

### Write cells dataframe
cell_df['idx'] = range(len(cell_df))
print(cell_df)
cell_df.to_pickle("cells.pkl")


d_x = 0.01

ax = plt.gca()

r_pixel  = [ 39, 75, 155, 213, 271 ]
dz_pixel  = 280

r_strips = [ 405, 562, 762, 1000 ]
dz_strips = 1150

z_pixel_ec  = [ 350, 420, 530, 670, 870, 1100, 1400, 2000, 2300, 2650 ]
z_strips_ec = [ 1300, 1600, 1900, 2250, 2650 ]

plt.grid(True)

for i_r, r in enumerate( r_pixel ):
    label = None
    if i_r == 0:
        label = "Tracker"
    plt.plot( [0.0, dz_pixel], [r, r], color = colors[9], label = label, lw = linewidth )

for i_r, r in enumerate( r_strips ):
    plt.plot( [0.0, dz_strips], [r, r], color = colors[9], lw = linewidth )

for z in z_pixel_ec[:7]:
    plt.plot( [ z, z ], [ 39, 271 ],  color = colors[9], lw = linewidth )
    
for z in z_pixel_ec[7:]:
    plt.plot( [ z, z ], [ 155, 271],  color = colors[9], lw = linewidth )

for z in z_strips_ec:
    plt.plot( [z, z], [405, 1000], color = colors[9], lw = linewidth )

r_iron  = [ 1166.4433333333333, 1332.8866666666666 ]
dz_iron = [ 2527.27, 2839.51 ]

label = "Iron"
for r, dz in zip( r_iron, dz_iron ):
    plt.plot( [ 0.0, dz ], [ r, r ], color = colors[8], lw = linewidth, label = label )
    label = None

plt.plot( [ 2839.51, 2839.51 ], [ 292.797, 1332.8866666666666 ], color = colors[8], lw = linewidth, label = label )
        
plt.xlim( 0.0, 6e3 )
plt.ylim( 0.0, 6e3 )

plt.xlabel('z [mm]')
plt.ylabel('y [mm]')
plt.gca().set_axisbelow(True)
plt.legend()
plt.gca().legend( fancybox = 'round', frameon = True, loc = "upper right", fontsize = 5.5 )
plt.tight_layout()

plt.show()


