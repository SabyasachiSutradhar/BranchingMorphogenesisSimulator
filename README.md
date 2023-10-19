# BranchingMorphogenesisSimulator
C code written by Sabyasachi Sutradhar, Yale university (sabyasachi.sutradhar@yale.edu) 2022
Algorithm:
This code simulates a branched tree structure in two dimensions. Initially, 2-4 randomly oriented branches are generated from the origin (0,0) and then the tips go through three state dynamics G-P-S with different rates Kgp (signifying transition rate from Growth to shrinkage phase), Kgs, Kpg, Kps, Ksg, Ksp as provided by the user. New branches spawned from an existing branch at a random place with a user-provided rate. For a more detailed description and algorithm please see Shree et al, Sci. Adv, 2022.  A file containing the following parameters is necessary to run the code:

///////Parameters for simulating branched neuronal Morphology
N_Sample=1;///Number of samples to simulate
Time_Steps=10000;////Maximum number of time steps, each time step is 0.1
Tip_Persis=250.0;////Persistence length of soma in um
BranchingAngleMean=90.0;////Mean branching angle wrt mother branch in degrees
BranchingAngleStd=5.0;////Std dev of branching angle wrt mother branch in degrees
BranchingRate=0.0025;////branching rate per min per micron;
Vg=1.0;///Tip Growth Velocity micron/min
Vs=1.0;//absolute value of tip shrinkage velocity micron/min
Kgp=0.50;//Transition rate from growth to pause /min
Kgs=0.25;//Transition rate from growth to shrinkage /min
Kpg=0.50;//Transition rate from pause to growth /min
Kps=0.50;//Transition rate from pause to shrinkage /min
Ksg=0.50;//Transition rate from shrinkage to pause /min
Ksp=0.50;//Transition rate from shrinkage to growth /min
//////////////// Boundary ////////////////////////////////
Boundary_type=repulsive;//options are "repulsive", "static", "free" and "pbc"
SizeX=200.0;//// Box  Size in X direction in micron
SizeY=200.0;//// Box Size in Y direction in micron
//////////////////////////////// Data dumping parameters /////////////////////////
Print_Conf=yes;///Print Configuration Images in pgm format (write yes/no)
Dump_Conf=1000;///Dump image every 
Dump_Data=1;///Dump image every 
pixelsize=0.1;////in micron

Output:
The code outputs neuron morphology data in .swc format along with skeletonized image in .pgm format. Additionally, the code produce a Timeseries data containing Branch number, branch lengths etc.

Requirements to run the code:
1.	gcc
2.	cmake (VERSION 3.24)
3.	UNIX/Mac Operating System

Building and installing the code:
To build the code run the following commands sequentially.

cd  “path to NeuronMorphologySimulator”
cmake -S src/ -B build/
cmake --build build/
cd build
make install  (you might need to use sudo make install if the access is denied)

Running the code
To run the code just type the following in the command window:
NMorphSim Parameters.in

Example of an output image file:

![ConfIm-Sample-1-TimeStep-10000](https://github.com/SabyasachiSutradhar/BranchingMorphogenesisSimulator/assets/49563656/66fec774-d052-48f9-a5d5-3513bd1b3d50)


