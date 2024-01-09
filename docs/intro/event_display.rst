Phoenix event display
---------------------

**Live demo hosted here:** https://cocoa-phoenix.web.cern.ch/

The ``phoenix`` directory contains the ingredients for displaying COCOA
events using the `HSF Phoenix
software <https://github.com/HSF/phoenix>`__.

-  ``event`` subdirectory: scripts for dumping COCOA output ROOT files
   into the suitable json format.
-  ``packages`` subdirectory: the changed files with respect to the
   Phoenix repository, with directory structure preserved. Note that 
   this builds the COCOA geometry.

Steps to get it fired up:

1. clone and follow the README on the Phoenix repository to get it set
   up locally.
2. replace the cloned files with the ones in the
   ``COCOA/phoenix/packages``. Note that this has only been tested at `a
   specific snapshot <https://github.com/HSF/phoenix/pull/536>`__ in the
   Phoenix code history.
3. (this step can be skipped in favor of using the default event files
   provided). Use the ``dump_phoenix_eventdata.py`` script to parse a
   COCOA output file, for example:

    .. code-block:: none
    
            python phoenix/event/dump_hdf5.py -i path/to/input_COCOA_file.root -o path/to/output_event_file.json -n 1

4. Copy the json event file to
   ``packages/phoenix-ng/projects/phoenix-app/src/assets/files/cocoa/``
   and edit the ``eventFile`` field in
   ``packages/phoenix-ng/projects/phoenix-app/src/app/sections/cocoa/cocoa.component.ts``
   appropriately.
5. Compile phoenix with yarn and open in browser window!
