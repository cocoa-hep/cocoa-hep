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
        if branch_name in ["cell_parent_list","particle_to_node_idx","particle_to_node_weight","cell_parent_energy","topo_jet_constituents_list","true_jet_constituents_list","supercluster_topos"]:
            continue
        arr = tree[branch_name].array(library="np", entry_stop=lastevent,entry_start=firstevent)
        data_array[branch_name] = arr
    
    return data_array

def pad_branch(branch_array):

    lengths = [event.shape[0] for event in branch_array]
    pad_with = np.NaN if isinstance(branch_array.flat[0], np.floating) else -9999
    output_list = []
    padded_branch_array = np.stack([
        np.pad(
            event, 
            (0, max(lengths) - lengths[idx]),
            "constant",
            constant_values = pad_with
        )
        for idx, event in enumerate(branch_array)
    ])

    return padded_branch_array

def writer(output_path,data_array,save_jagged=True):

    file = h5py.File(output_path, 'w')

    print("Writing branches...")
    for branch_name, branch_array in tqdm(data_array.items()):

        if not isinstance(branch_array[0],np.ndarray):
            continue

        num_entries = len(branch_array)

        if save_jagged: # see https://docs.h5py.org/en/stable/special.html#arbitrary-vlen-data
            dt = h5py.vlen_dtype(np.dtype('float32'))
            shape = (num_entries, )

            dataset = file.create_dataset(
                branch_name,
                shape,
                dtype=dt,
                compression="gzip",
                chunks=True)

            for idx,event in enumerate(branch_array):
                dataset[idx] = event
        else:
            dt = 'f'
            branch_array = pad_branch(branch_array)
            dataset = file.create_dataset(
                branch_name,
                branch_array.shape,
                data=branch_array,
                dtype=dt,
                compression="gzip",
                chunks=True)

    file.close()
    print("Output written to ",output_path)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",   dest="input",         type=str,  help="path to input ROOT file", required=True)
    parser.add_argument("-t","--tree",    dest="tree",          type=str,  help="name of TTree", default="")
    parser.add_argument("-o","--output",  dest="output",        type=str,  help="path to output h5 file", default="")
    parser.add_argument("-n","--nevents", dest="nevents",       type=int,  help="number of events to parse", default=-1)
    parser.add_argument("-s","--start",   dest="start",         type=int,  help="event to start on", default=0)
    parser.add_argument("-j","--jagged",  dest="save_jagged",   type=int,  help="save output as jagged array (alternative is padded array)", default=1)
    parser.set_defaults(save_jagged=False)
    args = parser.parse_args()

    if args.output == "":
        args.output = args.input.replace(".root",".h5")

    data_array = reader(args.input, args.nevents, args.start, args.tree)
    writer(args.output,data_array, bool(args.save_jagged))

if __name__ == '__main__':
    main()