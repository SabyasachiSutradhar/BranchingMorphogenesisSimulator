
/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header executes the dynamic memory aloocation of relevant variables.
Copyright @ Sabyasachi Sutradhar
*/

/////////////////////////////////
Tipx=dmatrix(0,N_TipMax,0,Max_Tippoints+1);
Tipy=dmatrix(0,N_TipMax,0,Max_Tippoints+1);
ParentID=imatrix(0,N_TipMax,0,Max_Tippoints+1);
tippoints=dmatrix(0,Max_Tippoints,0,2);
Tipl=ivector(0,N_TipMax);
Tiptype=ivector(0,N_TipMax);
tipx=dvector(0,Max_Tippoints);
tipy=dvector(0,Max_Tippoints);
Tiplength=dvector(0,N_TipMax);
topneighbour=ivector(0,N_TipMax);
bottomneighbour=ivector(0,N_TipMax);
topneighbourindex=imatrix(0,N_TipMax,0,2);
bottomneighbourindex=imatrix(0,N_TipMax,0,2);
new_tipindex=ivector(0,N_TipMax);
state=ivector(0,N_TipMax);
v=dvector(0,N_TipMax);
tipendx=dvector(0,N_TipMax);
tipendy=dvector(0,N_TipMax);
