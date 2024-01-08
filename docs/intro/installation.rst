Installation
======

The COCOA repository can be downloaded from git repository `cocoa-hep <https://github.com/cocoa-hep/cocoa-hep.git>`_ . 

Please follow the steps below to install COCOA.

Install from Docker
--------
The most convenient way to install COCOA is to use its `Docker <https://github.com/cocoa-hep/cocoa-hep/blob/main/Dockerfile>`_ image from the home area. 

After installing docker execute the following commands
    .. code-block:: none
    
            docker pull ghcr.io/cocoa-hep/cocoa-hep:main
            docker image tag $(docker images | grep cocoa-hep | head -n 1 | awk '{print $3}') cocoa-hep
            docker run -it cocoa-hep

Non-Docker installation
--------

To simplest way to prepare all dependencies is to mount the `CernVM File System <https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html>`_  and run
    .. code-block:: none
    
            source cocoa-hep/setup_cvmfs.sh

After the environment setup, the following commands in the shell must return the environments properly : 
    .. code-block:: none
    
            geant4-config --cflags
            root-config --cflags
            pythia8-config --cflags

Then in the `COCOA <https://github.com/cocoa-hep/cocoa-hep/tree/main/COCOA>`_ directory run the following commands:
    .. code-block:: none
    
            mkdir build
            cd build
            cmake ../
            make -j<# cpu cores>
            cd ..

Run
--------
From within `COCOA <https://github.com/cocoa-hep/cocoa-hep/tree/main/COCOA>`_ directory:

        .. code-block:: rst 

           **./build/COCOA** run with Geant4 User Interface.
           **./build/COCOA -h**  show input options for batch-mode with following options
           **--config (-c) <str>** path to json configuration file
           **--macro (-m) <str>** path to Geant4 or Pythia8 macro file for event generation (can be set in json configuration file)
           **--output (-o) <str>** path (incl. name) of output ROOT file to be written (can be set in json configuration file)
           **--input (-i) <str>** path to HepMC (.hmc) input file (overrides the default path set in the HepMC macro file)
           **--seed (-s) <int>** set random seed
           **--nevents (-n) <int>** number of events to generate (default is read from macro)
Example
-------- 
An example to run the code interactively:

        .. code-block:: none 

           ./build/COCOA --macro  /path/to/COCOA/COCOA/macro/Pythia8/ttbar.in --config  /path/to/COCOA/COCOA/config/config_doc.json  /path/to/outputdir/output_name.root --seed 5

Convert
-------- 
To convert the output files from COCOA from ROOT to hdf5 format, the `util/dump_hdf5.py` can be used as follows:

        .. code-block:: none 

            python util/dump_hdf5.py -i path/to/input.root -o path/to/output.h5

To see more options, pass the `-h` argument.
