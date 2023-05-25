<img src="https://github.com/scd-hep/scd-hep/blob/main/docs/imgs/cocoa.png" height="250">

# COCOA
The COnfigurable Calorimeter simulatiOn for Ai (COCOA) is a nearly-hermetic calorimeter simulated with [Geant4](https://geant4.web.cern.ch) and interfaced to the [Pythia8](https://www.pythia.org) event generator. This open-source simulation is aimed to support the development of machine learning algorithms in high energy physics that rely on realistic particle shower modeling, such as reconstruction, fast simulation, and low-level analysis.

The COCOA calorimeter comprises a barrel and endcap system with configurable granularity, and with nearly uniform material depth distribution in pseudorapidity. An inner tracker system consisting of silicon and iron layers immersed in a magnetic field can be included optionally, along with basic tracking emulation. Output data are processed using on-board algorithms for topological clustering of calorimeter cells, graph creation, and jet clustering. The COCOA geometry is also provided in a format supporting event visualization with [Phoenix](https://github.com/HSF/phoenix).

## Publication
[![DOI](https://zenodo.org/badge/563008933.svg)](https://zenodo.org/badge/latestdoi/563008933)

[![Docs](https://img.shields.io/badge/Docs-cocoa--hep-green.svg)](https://cocoa-hep.readthedocs.io/en/latest/index.html) [![arXiv](https://img.shields.io/badge/arXiv-2303.02101-red.svg)](https://arxiv.org/abs/2303.02101)

## Install

### Docker

The most convenient way to install COCOA is to use its docker image:
```
docker pull ghcr.io/cocoa-hep/cocoa-hep:main
docker image tag $(docker images | grep cocoa-hep | head -n 1 | awk '{print $3}') cocoa-hep
docker run -it cocoa-hep
```

Please note that in this container cocoa and its dependencies are installed in /root .

### Non-Docker

To simplest way to prepare all dependencies is to mount the [CernVM File System](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html) and run
```
source setup_cvmfs.sh
```
Otherwise the dependencies need to be taken care of individually. The [Dockerfile](Dockerfile) can be used for guidance in this case.


Then in the `COCOA` directory run the following commands:
```
mkdir build
cd build
cmake ../
make -j<# cpu cores>
cd ..
```

## Run
From within `COCOA` directory:

`./build/COCOA` - run with Geant4 User Interface.

`./build/COCOA -h` - show input options for batch-mode.

**List of options:**
- `--config (-c) <str>` – path to json configuration file.
- `--macro (-m) <str>` – path to Geant4 or Pythia8 macro file for event generation (can be set in json configuration file.
- `--output (-o) <str>`  – path (incl. name) of output ROOT file to be written (can be set in json configuration file).
- `--seed (-s) <int>` –   set random seed.
- `--nevents (-n) <int>` - number of events to generate (default is taken from macro).


**Example:**
```
./build/COCOA --macro  /path/to/COCOA/COCOA/macro/Pythia8/ttbar.in --config  /path/to/COCOA/COCOA/config/config_doc.json  /path/to/outputdir/output_name.root --seed 5
```

## Convert
To convert the output files from COCOA from ROOT to hdf5 format, the `util/dump_hdf5.py` can be used as follows:
```
python util/dump_hdf5.py -i path/to/input.root -o path/to/output.h5
```
To see more options, pass the `-h` argument.

## Phoenix event display

<img src="https://github.com/scd-hep/scd-hep/blob/main/docs/imgs/ttbar.png" height="250">

The `phoenix` directory contains the ingredients for displaying COCOA events using the [HSF Phoenix software](https://github.com/HSF/phoenix). 
- `event` subdirectory: scripts for dumping COCOA output ROOT files into the suitable json format.
- `packages` subdirectory: the changed files with respect to the Phoenix repository, with directory structure preserved. Note that this builds the COCOA geometry.

Steps to get it fired up:
1. clone and follow the README on the Phoenix repository to get it set up locally.
2. replace the cloned files with the ones in the `COCOA/phoenix/packages`. Note that this has only been tested at [a specific snapshot](https://github.com/HSF/phoenix/pull/536) in the Phoenix code history.
3. (this step can be skipped in favor of using the default event files provided). Use the `dump_phoenix_eventdata.py` script to parse a COCOA output file, for example:
```
python phoenix/event/dump_hdf5.py -i path/to/input_COCOA_file.root -o path/to/output_event_file.json -n 1
```
4. Copy the json event file to `packages/phoenix-ng/projects/phoenix-app/src/assets/files/cocoa/` and edit the `eventFile` field in `packages/phoenix-ng/projects/phoenix-app/src/app/sections/cocoa/cocoa.component.ts` appropriately.
5. Compile phoenix with yarn and open in browser window!
