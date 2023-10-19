/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header file contains all the parameters necessary for simulating the branched dendritic arbor.
Copyright @ Sabyasachi Sutradhar
*/

void Time_IndependentParameters(){
////////////////////////////// Simulation parameters ////////////////////////////
Initial_TipLength=10.0;///initial tip length from soma in um
Tip_intv=0.1; ///tip increment length (lattice size) in um
Max_Tippoints=(int)(1000.0/Tip_intv);
Min_tiplength=0.50;//Minimum tip initiation length in um
BranchingLength_Threshold=0.80;//threshold length for branching
R_Branch=0.15;////average radius of the branches
Interaction_length=0.15;///collision interaction length between two tips in um
N_TipMax=100000;//Maximum number of tips
Dt=0.1; ///time step in min
R_Soma=5.0;////radius of Soma in um
///// Now include the user defined parametes
}
