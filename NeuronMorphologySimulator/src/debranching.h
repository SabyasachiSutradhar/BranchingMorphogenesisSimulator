/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header excutesdeletion of a branch and house keep the other relevant variables.
Copyright @ Sabyasachi Sutradhar
*/

int check_dynamic(int i,int j){
  int ind=1;
  if(Tiptype[i]==0 && Tiptype[j]==0 ){
    ind=0;
  }
  return ind;
}

////////  merge and assign top and bottom neighbours properly after branch deletion
void merge_tips(int i){
  double dis1,dis2,dist,dis;
  int indd,indt,mm,top1,top2,nd1,nd2,nt1,nt2;
  int ind1=bottomneighbourindex[i][1];
  int ind2=bottomneighbourindex[i][2];
  if(determine_downbranch(i,ind1)==1){indd=ind1;indt=ind2;}
  else{indd=ind2;indt=ind1;}
  mm=Tipl[indd];
  for(int j=2;j<=Tipl[indt]; j++){
    mm++;
    Tipx[indd][mm]=Tipx[indt][j];
    Tipy[indd][mm]=Tipy[indt][j];
  }
  Tipl[indd]=mm;
  Tiplength[indd]=calculate_curvelength(indd,Tipl[indd]);
  Tiptype[indd]=check_dynamic(ind1,ind2);
  if(Tiptype[indd]==0){
    topneighbour[indd]=1;
    top1=topneighbourindex[indt][1];top2=topneighbourindex[indt][2];
    allocate_top_neighbour(indd,top1,top2);
    allocate_bottom_neighbour(top1,top2,indd);
    allocate_bottom_neighbour(top2,top1,indd);
  }else{
    state[indd]=state[indt];
    v[indd]=v[indt];
    topneighbour[indd]=0;
  }
////////// delete one of the branches
  Tipl[indt]=0;
  Tiptype[indt]=0;
}
//////// debranching if the tip length is smaller than a threshold value
void debranching(){
  for (int i=1; i<=n_tip; i++){
    if(i>N_TipInitial){
      if(Tipl[i]==0 && Tiptype[i]==1){
        merge_tips(i);///merging two tip after a deletion event
      }
    }
  }
}
////////// relabel the branch indices after branch deletion
void relabel_index(){
  int ix,iy;
  int mm=0;
  ////// First find out the indices after branch removal
  for (int i=1; i<=n_tip; i++){
    if(Tipl[i]>0){
      mm++;
      new_tipindex[i]=mm;
    }
  }

n_tipend=0;
Total_Length=0.0;
//N_points=0.0;
mm=0;
  minx=1000000;maxx=-10000000;miny=10000000;maxy=-1000000;
  int k=0;
  /////////////////////////////////// Assign the new indices
  for (int i=1; i<=n_tip; i++){
    if(Tipl[i]>0){
      mm++;
      Tipl[mm]=Tipl[i];
      for (int j=1; j<=Tipl[mm]; j++){
        Tipx[mm][j]=Tipx[i][j];
        Tipy[mm][j]=Tipy[i][j];
/////////////// some info For configuration file ///////////////////
        ix=ceil(Tipx[mm][j]/pixelsize);iy=ceil(Tipy[mm][j]/pixelsize);
        if(ix<minx){minx=ix;}
        if(iy<miny){miny=iy;}
        if(ix>maxx){maxx=ix;}
        if(iy>maxy){maxy=iy;}
        if(j>=2){
          k+=1;
          ParentID[mm][j]=k;
        }
      }
///////////////////////////////////////// assign other tip properties
      Tiptype[mm]=Tiptype[i];
      Tiplength[mm]=Tiplength[i];
      Total_Length+=Tiplength[mm];
      topneighbour[mm]=topneighbour[i];
      bottomneighbour[mm]=bottomneighbour[i];
      if(Tiptype[mm]==1){
        n_tipend++;
        state[mm]=state[i];
        v[mm]=v[i];
        topneighbour[mm]=0;
        //////////////////// Hull Area_Hull  ////////////////////////////
        tipendx[n_tipend]=Tipx[mm][Tipl[mm]];
        tipendy[n_tipend]=Tipy[mm][Tipl[mm]];
      }
      if(topneighbour[mm]==1){
        topneighbourindex[mm][1]=new_tipindex[topneighbourindex[i][1]];
        topneighbourindex[mm][2]=new_tipindex[topneighbourindex[i][2]];
      }
      if(bottomneighbour[mm]==1){
        bottomneighbourindex[mm][1]=new_tipindex[bottomneighbourindex[i][1]];
        bottomneighbourindex[mm][2]=new_tipindex[bottomneighbourindex[i][2]];
      }
    }
  }

  n_tip=mm;
  ////////////////  printf("###############\n"); //////////////////////
}
