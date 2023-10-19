/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header is collection of necessary functions.
Copyright @ Sabyasachi Sutradhar
*/
double assign_direction(int i,int mm){
  double dir=wrap_angle(atan2(Tipy[i][mm]-Tipy[i][mm-1],Tipx[i][mm]-Tipx[i][mm-1]));
  return dir;
}

double calculate_curvelength(int i,int N){
  double dis=0.0;
  if(N>=2){
  for (int j=1; j<=N-1; j++){
    dis+=distance2D(Tipx[i][j],Tipy[i][j],Tipx[i][j+1],Tipy[i][j+1]);
  }
}else{
  dis=0.0;
}
  return dis;
}


void allocate_bottom_neighbour(int ind1,int ind2, int ind3){
    if(bottomneighbour[ind1]==1){
      int indmin=find_min(ind2,ind3);
      int indmax=find_max(ind2,ind3);
    bottomneighbourindex[ind1][1]=indmin;bottomneighbourindex[ind1][2]=indmax;
  }
}


void allocate_top_neighbour(int ind1,int ind2, int ind3){
  if(topneighbour[ind1]==1){
    int indmin=find_min(ind2,ind3);
    int indmax=find_max(ind2,ind3);
    topneighbourindex[ind1][1]=indmin;topneighbourindex[ind1][2]=indmax;
  }
}


int determine_downbranch(int i,int n1){
  int ind=0;
  if(topneighbour[n1]==1){
  if(topneighbourindex[n1][1]==i || topneighbourindex[n1][2]==i ){ind=1;}
}
  return ind;
}
