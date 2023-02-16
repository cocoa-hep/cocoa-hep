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
source SCD/setup.sh
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
