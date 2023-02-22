from tqdm import tqdm
import argparse
import uproot
import numpy as np
import h5py

def reader(input_path, nevents=-1, firstevent=0, tree_name=""):
    
    rootfile = uproot.open(input_path)
    if len(rootfile.keys()) == 0:
        raise Exception(f"Input file {input_path} contains no keys")
    if tree_name=="":
        tree_name = rootfile.keys()[0]
    print(f"Reading tree {tree_name} ...")
    tree = rootfile[tree_name]

    if nevents<0:
        nevents = tree.num_entries
    lastevent = min(tree.num_entries, firstevent + nevents)

    data_array = {}

    for branch_name in tree.keys():
        arr = tree[branch_name].array(library='np', entry_stop=lastevent,entry_start=firstevent)
        if isinstance(arr[0],np.ndarray):
            data_array[branch_name] = arr
    
    return data_array

def writer(output_path,data_array,save_jagged=True):

    file = h5py.File(output_path, 'w')

    if(save_jagged): # see https://docs.h5py.org/en/stable/special.html#arbitrary-vlen-data
        for branch_name, branch_array in tqdm(data_array.items()):

            #if branch_name in ["cell_parent_list","particle_to_node_idx","particle_to_node_weight","cell_parent_energy","topo_jet_constituents_list","true_jet_constituents_list"]:
            #    continue
            
            num_entries = len(branch_array)
            dt = h5py.vlen_dtype(np.dtype('float32'))
            ex = branch_array[0]
            shape = (num_entries,branch_array[0].shape[0]) #if len(branch_array[0].shape) == 1 else (num_entries,branch_array[0].shape[0],branch_array[0].shape[1])
            dataset = file.create_dataset(branch_name,shape,dtype=dt)
            for idx,row in enumerate(branch_array):
                dataset[idx] = row

    file.close()
    print("Output written to ",output_path)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",   dest="input",   type=str, help="path to input ROOT file", required=True)
    parser.add_argument("-t","--tree",    dest="tree",    type=str, help="name of TTree", default="")
    parser.add_argument("-o","--output",  dest="output",  type=str, help="path to output json file", default="events.json")
    parser.add_argument("-n","--nevents", dest="nevents", type=int, help="number of events to parse", default=-1)
    parser.add_argument("-s","--start",   dest="start",   type=int, help="event to start on", default=0)
    args = parser.parse_args()

    data_array = reader(args.input, args.nevents, args.start, args.tree)
    writer(args.output,data_array)

if __name__ == '__main__':
    main()