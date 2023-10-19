/*
C code written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
Copyright @ Sabyasachi Sutradhar
*/
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <stdarg.h>
  #include <stdbool.h>
  #include <stddef.h>
  
  static long idum=-123456797;
  #include "ran_generator.h"
  double pi=3.14159265;
  double PI=3.14159265;
  #include "nrutils.h"
  #define NUMELEM(a) (sizeof(a) / sizeof(*a))
/// @USER defined parameters ////////////////////////////////////////////////////
int N_Sample,Dump_Conf,Time_Steps,Dump_Data;
double Tip_Persis,BranchingAngleMean,BranchingAngleStd,BranchingRate,Vg,Vs,Kgp,Kgs,Kpg,Kps,Ksg,Ksp,SizeY,SizeX,pixelsize;
char Boundary_type[100],Print_Conf[100];
/// @Simulation internal parameters /////////////////////////////////////////////////
double time,Tip_intv,Dt,Total_Length,R_Soma,Initial_TipLength,R_Branch,BranchingLength_Threshold,MinBrLength,Min_tiplength,Interaction_length;
int N_TipMax,N_TipInitial,Max_Tippoints;
double **Tipx,**Tipy,**tippoints,*tipx,*tipy,*Tiplength;
int *state,*Tipl,*Tiptype,*topneighbour,*bottomneighbour,**topneighbourindex,**bottomneighbourindex,*new_tipindex,**ParentID;
int n_tip,n_tipend,t_max,n_pts,new_tip;
double *tipendx,*tipendy,*v;
//////////////////////////////////////////
int n_im,mod_data,image_counter,tint;
int minx,maxx,miny,maxy;
//////////////////////// File Configuration ////////////////////////////
  char fileconf[1000],file1[1000],fileconfim[1000];
  FILE *fpconf,*fpconfim,*fp1;
  /////////////////////////////////////////////////////
////////////////////////////// essential subroutines
  #include "ReadInput.h"
  #include "Simulation_Parameters.h"
  #include "essentials.h"
  #include "ContactBasedRetraction.h"
  #include "debranching.h"
  #include "tip_dynamics.h"
  #include "branch_formation.h"
  #include "Initialization.h"
  #include "Print_Data.h"
//////////////////////////////////////////////////
/////////////////// Main function //////////////////
////////////////////////////////////////////////////
  int main(int argc, char *argv[]){
    ////////////////// Read parameters from file
int tt;
tt=ReadInputData(argc,argv[0],argv[1]);
    if(tt == -1) {
		perror("Error");
		exit(1);
	}
  ///////////////////////////////
    double time_h;
  /////////////////// Simulation and static parameters /////////
  Time_IndependentParameters();
  ///////////////////////////////////////////////////////////////////
  #include "mem_allocation.h"//allocate memory for variables

    /////////////////// SAMPLE LOOP //////////////////////////
    for(int sample=1;sample<=N_Sample;sample++){
      sprintf(file1,"Time_NumberofTips-Sample-%d.csv",sample);//filename for each sample
      fp1=fopen(file1,"w");
      fprintf(fp1,"#Simulated data for input file %s\n",argv[1]);
      fprintf(fp1,"Timestep, Number of branches, Number of tips,Total_Length (micron)\n");
      fclose(fp1);
  /////////////////////// Initialization /////////////////////
      initialize();
  /////////////////////// Initial time  ///////////////////////////
time=0.0;
/////////////////////// Time loop stats here/////////////////////////////////
for(tint=0;tint<=Time_Steps;tint++){
  //////////////////////////////////////////////////////////////////////////////
            Branching();
            Tip_Dynamics();
            debranching();
            relabel_index();
  //////////////// Dump Data and Confs//////////////////////////////////////////////////
      if(strcmp(Print_Conf,"yes")==0 && tint%Dump_Conf==0 ){print_configurations(sample);}
      if(tint%Dump_Data==0){print_datafiles();}
  ///////////////////////////////////////////////////////////////////////////////////
    time=time+Dt;
    }////////////////////// End of time loop ///////////////////////
  } /////////////////////end of sample loop//////////////////
  return 0;
}
 




