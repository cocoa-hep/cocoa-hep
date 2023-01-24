Installation
======

Please follow the steps below to install COSA.

Install from Docker
--------
The most convenient way to install COSA is to use its `Docker <https://github.com/scd-hep/scd-hep/blob/main/Dockerfile>`_ image from the home area. 

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
           **-path_to_config** path to json configuration file
           **-path_to_script** path to Geant4 macro file for particle gun (can be set in json configuration file)
           **-path_to_output** destination and name of the output root file (can be set in json configuration file)
           **-set_seed_value** set random seed
Example
-------- 
An example to run the code interactively:

        .. code-block:: none 

           ./build/SCDMain -path_to_script <path_to_COSA>/SCD/macro/Pythia8/ttbar.in   -path_to_config <path_to_COSA>/SCD/config/config_doc.json   -path_to_output <path_to_output>output_name.root -set_seed_value 1

