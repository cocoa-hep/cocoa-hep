Installation
======

Please follow the steps below to install COSA.

Install from Docker
--------
The most convenient way to install SCD is to use its `Docker <https://github.com/scd-hep/scd-hep/blob/main/Dockerfile>`_ image from the home area. 

After installing docker execute the following commands
    .. code-block:: none
    
            docker pull ghcr.io/scd-hep/scd-hep:main
            docker image tag $(docker images | grep scd-hep | head -n 1 | awk '{print $3}') scd-hep
            docker run -it scd-hep

Non-Docker installation
--------
The COSA repository can be downloaded from git `scd-hep <https://github.com/scd-hep/scd-hep.git>`_

To simplest way to prepare all dependencies is to mount the `CernVM File System <https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html>`_  and run
    .. code-block:: none
    
            source scd-hep/setup_cvmfs.sh

After the environment setup, the following commands in the shell must return the environments properly : 
    .. code-block:: none
    
            geant4-config --cflags
            root-config --cflags
            pythia8-config --cflags

Then in the `SCD <https://github.com/scd-hep/scd-hep/tree/main/SCD>`_ directory run the following commands (see `setup_cvfms.sh <https://github.com/scd-hep/scd-hep/blob/main/setup_cvmfs.sh>`_ ):
    .. code-block:: none
    
            mkdir build
            cd build
            cmake ../
            make -j<# cpu cores>
            cd ..

Run
--------
From within `SCD <https://github.com/scd-hep/scd-hep/tree/main/SCD>`_ directory:

        .. code-block:: rst 

           **./build/SCDMain** run with Geant4 User Interface.
           **./build/SCDMain -h**  show input options for batch-mode with following options
           **--config (-c) <str>** path to json configuration file
           **--macro (-m) <str>** path to Geant4 or Pythia8 macro file for event generation (can be set in json configuration file)
           **--output (-o) <str>** path (incl. name) of output ROOT file to be written (can be set in json configuration file)
           **--seed (-s) <int>** set random seed
           **--nevents (-n) <int>** number of events to generate (default is read from macro)
Example
-------- 
An example to run the code interactively:

        .. code-block:: none 

           ./build/SCDMain --macro  /path/to/SCD/SCD/macro/Pythia8/ttbar.in --config  /path/to/SCD/SCD/config/config_doc.json  /path/to/outputdir/output_name.root --seed 5

Convert
-------- 
To convert the output files from SCD from ROOT to hdf5 format, the `util/dump_hdf5.py` can be used as follows:

        .. code-block:: none 

            python util/dump_hdf5.py -i path/to/input.root -o path/to/output.h5

To see more options, pass the `-h` argument.