Gerasimova-Zatsepin effect simulation
=====================================

Code by Lafèbre for GZ effect simulations


How to run the GZ Simulation
----------------------------

1. `magneticfield_1.cpp` is a file that can be run independently from
the other files. It can be used to plot the different components of the
interplanetary magnetic field.

2. In order to run the other files, they have to be put in one folder.
All the .cpp files have to be compiled into one .out file. It is not
possible to compile `main_.cpp` and `id_.cpp` while they are in the same
folder.

3. The `main_.cpp` can be used to create a database of GZ-events, which
are deflected by the interplanetary magnetic field.

4. The database is written to a .txt file which can be used as input for
`id_.cpp`.

5. The `id_.cpp` can then be used to analyse the database of particles
and implement a detector array.

6. The detector specifics are defined in `Detector.cpp` and
`DetHisparc.cpp`.

7. The .dat files that are created in `id_.cpp` can be interpreted using
the information from `Hist2D.cpp` and `id_.cpp`.
