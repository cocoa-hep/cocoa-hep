from tqdm import tqdm
import argparse
import uproot
import numpy as np
import h5py

def reader(input_path, nevents=-1, firstevent=0, tree_name="", light=False):
    
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

        arr = tree[branch_name].array(library="np", entry_stop=lastevent,entry_start=firstevent)

        if isinstance(arr[0],uproot.STLVector):

            if light: #skip branches containing 2D arrays
                continue

            arr = np.array([ #branch
                        np.array([ #event
                            np.array( #row
                                np.array(row.tolist() if len(row)!=0 else [np.NaN]) #column
                            )
                            for row in event
                        ],dtype=object)
                    for event in arr
                ],dtype=object)

        data_array[branch_name] = arr
    
    return data_array

def pad_array(arr,_max_length = (-1,-1)):

    if arr[0].dtype==object and hasattr(arr[0][0], '__iter__'): #array of 2-dimensional arrays
        pad_with = np.NaN if arr[0][0].dtype.kind == 'f' else -9999
        lengths_0 = [len(event) for event in arr]
        lengths_1 = [ [len(node) for node in event] for event in arr]
        max_length_0 = max(lengths_0) if _max_length[0] < 0 else _max_length[0]
        max_length_1 = max([max(node_lengths) if len(node_lengths)!=0 else 0 for node_lengths in lengths_1]) if _max_length[1] < 0 else _max_length[1]

        padded_arr = np.stack([
            np.pad(
                pad_array(event,(max_length_1,-1)),
                [(0, max_length_0 - lengths_0[idx]),(0,0)],
                "constant",
                constant_values = pad_with
            )
            for idx, event in enumerate(arr)
        ]).astype('float32')

    else:
        pad_with = np.NaN if arr[0].dtype.kind == 'f' else -9999
        lengths = [len(event) for event in arr]
        max_length = max(lengths) if _max_length[0] < 0 else _max_length[0]

        padded_arr = np.stack([
            np.pad(
                event, 
                (0, max_length - lengths[idx]),
                "constant",
                constant_values = pad_with
            )
            for idx, event in enumerate(arr)
        ])

    return padded_arr

def writer(output_path,data_array,save_jagged=True):

    file = h5py.File(output_path, 'w')

    print("Writing branches...")
    for branch_name, branch_array in tqdm(data_array.items()):
        num_entries = len(branch_array)

        # for saving jagged arrays, see https://docs.h5py.org/en/stable/special.html#arbitrary-vlen-data
        # (not supported for branches storing 2D arrays)
        if save_jagged and branch_array[0].dtype!=object:
            shape = (num_entries, )
            dt = h5py.vlen_dtype(np.dtype('float32'))

        else:
            branch_array = pad_array(branch_array)
            shape = branch_array.shape
            dt = type(branch_array[0][0][0])

        dataset = file.create_dataset(
            branch_name,
            shape,
            dtype=dt,
            compression="gzip",
            chunks=True)

        for idx,event in enumerate(branch_array):
            dataset[idx] = event

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
    parser.add_argument("-l","--light",   dest="light",         type=int,  help="do not save expensive branches containing 2D arrays", default=0)
    args = parser.parse_args()

    if args.output == "":
        args.output = args.input.replace(".root",".h5")

    data_array = reader(args.input, args.nevents, args.start, args.tree, bool(args.light))
    writer(args.output,data_array, bool(args.save_jagged))

if __name__ == '__main__':
    main()