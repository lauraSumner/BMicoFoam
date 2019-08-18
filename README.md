# BMicoFoam
OpenFOAM 3.0.1 solver to model the fluid mechanics of the cochlea.
In order to use this code, you need to define a mesh with an internal wall (this will be the Basilar Membrane). This is possible through the createBaffles command. You will need to create a dictionary file in order to do this! Which you need to save in the case directory. 
The code is written to include 5000 particles in order to track their motion. You must therefore also create a particleProperites dictionary which contains the initial coordinates of each of the 5000 particles. 
There are also 2 more walls which move, called the oval window and the round window. The solver reads the coordinates 
of the BM, oval and round windows on each time step, so these must be calculated beforehand (see the WKB_travellingWave 
section for a MATlAB code to do that). The filenames for these coordinates must be as follows:

BM Deflection Coordinates and dimensions files:
matrix_x.txt matrix_y.text dimension_y.txt

Oval Window Coordinates:
ovalXcoords.txt ovalYcoords.txt

Round Window Coordinates:
roundXcoords.txt roundYcoords.txt

You must also label the BM membrane BM1.

createBafflesDict and particleProperties files in the main case file you should be able to run the code as long as the mesh is correct.

