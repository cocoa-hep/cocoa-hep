# Simplified Cylindrical Detector

Simplified Cylindrical Detector (SCD) is Geant4 and Pythia8 based simulation of the calorimetric system in the barrel and endcap region of ATLAS-like detector with postprocessing 
algorithms like topological clustering, particle flow and jet clustering algorithms. 

For more details, see [wiki](https://gitlab.com/anton70406/master/-/wikis/Simplified-Cylindrical-Detector).

**Requirements**:
- Geant4 >= 10.05
- ROOT >= 6.18
- HepMC >= 2.06
- jsoncpp >= 1.8
- cmake >=3.14
- Pythia8 >=8301
- FASTJet >= 3.3.2
## Install
### Download
**N.B.** the repo has migrated from AntonCh-G to wisroma-pflow.

Clone the repository:
`git clone git@github.com:wisroma-pflow/SCD.git`

### How To Build

Inside `SCD` directory: (See `make.sh`)
```
mkdir build
cd build
cmake ../
make
cd ..
```

## First run
From within `SCD` directory:

`./build/SCDMain` - run with Geant4 User Interface.

`./build/SCDMain -h` - show input options for batch-mode.

**List of options:**
- `-path_to_config <path>` – path to json configuration file.
- `-path_to_script <path>` – path to Geant4 macro file for particle gun (can be set in json configuration file).
- `-path_to_output <path>`  – destination and name of the output root file (can be set in json configuration file).
- `-set_seed_value <int>` –   set random seed.

**Example:**
```
./build/SCDMain -path_to_script  /path/to/SCD/SCD/macro/Pythia8/ttbar.in -path_to_config  /path/to/SCD/SCD/config/config_doc.json  /path/to/outputdir/output_name.root 
-set_seed_value 5
```
