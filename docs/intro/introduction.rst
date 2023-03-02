Introduction
======

The CalOrimeter Simulation for Ai applications (COSA) is a simulated calorimeter system based on Geant4. 
Primary particles can be generated with Pythia8 within COSA. 
The detector consists of a barrel and endcap calorimeter system with adjustable granularity. 
An inner detector consisting of silicon and iron layers can be included as an option. 
The detector layout is similar to the one of the ATLAS detector. 
Postprocessing algorithms, namely topological clustering of calorimeter cells, 
particle flow and jet clustering algorithms are provided.

This open-source simulation is aimed
to support the development of machine learning algorithms
in high energy physics that rely on realistic particle shower
descriptions, such as reconstruction, fast simulation, and low-level analysis. 
The package also includes Phoenix event display for visualizing events.

.. image:: ../imgs/ttbar.png
   :width: 100%
