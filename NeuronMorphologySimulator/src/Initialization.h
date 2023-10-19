/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header initialize the branches and the segment in which it grows.
Copyright @ Sabyasachi Sutradhar
*/

void initialize(){
///////////////////// Initialization //////////////////////
N_TipInitial=random_int(2,5); //number of initial branches
n_tip=N_TipInitial;
n_tipend=n_tip;
for (int i=1; i<=N_TipInitial; i++){
  Tiptype[i]=1;
  topneighbour[i]=0;
  bottomneighbour[i]=0;
  state[i]=1;
  v[i]=Vg;
  double dir=wrap_angle(gaussdev(0.0,pi/12.0)+(double)i*2.0*pi/(N_TipInitial));
  tipx[1]=0.0;tipy[1]=0.0;
  tipx[2]=Initial_TipLength*cos(dir);tipy[2]=Initial_TipLength*sin(dir);
  Tiplength[i]=Initial_TipLength;
  n_pts=(int) ceil(Tiplength[i]/Tip_intv);
  Tipl[i]=n_pts;
  tippoints=equispace_points(tipx,tipy,2,n_pts);
  for(int j=1;j<=n_pts;j++){
    Tipx[i][j]=tippoints[j][1];Tipy[i][j]=tippoints[j][2];
  }
}

}
