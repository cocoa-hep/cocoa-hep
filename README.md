# scd-hep
SCD Super Cool Detector

The Simplified Cylindrical Detector (SCD) is a simulated calorimeter system based on Geant4. Primary particles can be generated with Pythia8 within SCD. The detector consists of a barrel and endcap calorimeter system with adjustable granularity. An inner detector consisting of silicon and iron layers can be included as an option. The detector layout is similar to the one of the ATLAS detector. Postprocessing algorithms, namely topological clustering of calorimeter cells, particle flow and jet clustering algorithms are provided.

For more details, see [wiki](https://gitlab.com/anton70406/master/-/wikis/Simplified-Cylindrical-Detector).

## Install

### Docker

The most convenient way to install SCD is to use its docker image from <TO BE PROVIDED ONCE THE IMAGE IS PUBLIC>.

### Non-Docker

To simplest way to prepare all dependencies is to mount the [CernVM File System](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html) and run
```
source setup_cvmfs.sh
```
Otherwise the dependencies need to be taken care of individually. The [Dockerfile](Dockerfile) can be used for guidance in this case.


Then in the `SCD` directory run the following commands (see `make.sh`):
```
mkdir build
cd build
cmake ../
make -j<# cpu cores>
cd ..
```

## Run
From within `SCD` directory:

`./build/SCDMain` - run with Geant4 User Interface.

`./build/SCDMain -h` - show input options for batch-mode.

**List of options:**
- `--config (-c) <str>` – path to json configuration file.
- `--macro (-m) <str>` – path to Geant4 or Pythia8 macro file for event generation (can be set in json configuration file.
- `--output (-o) <str>`  – path (incl. name) of output ROOT file to be written (can be set in json configuration file).
- `--seed (-s) <int>` –   set random seed.
- `--nevents (-n) <int>` - number of events to generate (default is taken from macro).


**Example:**
```
./build/SCDMain --macro  /path/to/SCD/SCD/macro/Pythia8/ttbar.in --config  /path/to/SCD/SCD/config/config_doc.json  /path/to/outputdir/output_name.root --seed 5
```

## Convert
To convert the output files from SCD from ROOT to hdf5 format, the `util/dump_hdf5.py` can be used as follows:
```
python util/dump_hdf5.py -i path/to/input.root -o path/to/output.h5
```
To see more options, pass the `-h` argument.

## Phoenix event display

The `phoenix` directory contains the ingredients for displaying SCD events using the [HSF Phoenix software](https://github.com/HSF/phoenix). 
- `event` subdirectory: scripts for dumping SCD output ROOT files into the suitable json format.
- `packages` subdirectory: the changed files with respect to the Phoenix repository, with directory structure preserved. Note that this builds the SCD geometry.

Steps to get it fired up:
1. clone and follow the README on the Phoenix repository to get it set up locally.
2. replace the cloned files with the ones in the `SCD/phoenix/packages`. Note that this has only been tested at [a specific snapshot](https://github.com/HSF/phoenix/pull/536) in the Phoenix code history.
3. (this step can be skipped in favor of using the default event files provided). Use the `dump_phoenix_eventdata.py` script to parse a SCD output file, for example:
```
python phoenix/event/dump_hdf5.py -i path/to/input_SCD_file.root -o path/to/output_event_file.json -n 1
```
4. Copy the json event file to `packages/phoenix-ng/projects/phoenix-app/src/assets/files/scd/` and edit the `eventFile` field in `packages/phoenix-ng/projects/phoenix-app/src/app/sections/scd/scd.component.ts` appropriately.
5. Compile phoenix with yarn and open in browser window!