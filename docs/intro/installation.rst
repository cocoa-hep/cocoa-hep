Installation
======

Please follow the steps below to install SCD-HEP.

Install from Docker
--------
The most convenient way to install SCD is to use its `Docker <https://github.com/scd-hep/scd-hep/blob/main/Dockerfile>`_ image from the home area. 

Non-Docker installation
--------

To simplest way to prepare all dependencies is to mount the `CernVM File System <https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html>`_  and run
    .. code-block:: none
    
            source SCD/setup.sh
            
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

        .. code-block:: none 

           ./build/SCDMain - run with Geant4 User Interface.
           ./build/SCDMain -h - show input options for batch-mode.
           
Example
-------- 
An example to run the code interactively:

        .. code-block:: none 

           ./build/SCDMain -path_to_script  /path/to/SCD/SCD/macro/Pythia8/ttbar.in -path_to_config  /path/to/SCD/SCD/config/config_doc.json  /path/to/outputdir/output_name.root -set_seed_value 5

