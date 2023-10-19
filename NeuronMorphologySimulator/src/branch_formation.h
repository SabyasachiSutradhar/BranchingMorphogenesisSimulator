/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header puts a new branch randomly generated from an existing branch.
Copyright @ Sabyasachi Sutradhar
*/
double sample_branching_angle(double Br_dir){
///// Truncate branching wrap_angle
  double theta_a=pi/90.0;
double th=gaussdev(BranchingAngleMean,BranchingAngleStd);
  while(th<theta_a || th>pi-theta_a){
    th=gaussdev(BranchingAngleMean,BranchingAngleStd);
  }
  double thf=wrap_angle(Br_dir+random_sign()*th);
  return(thf);
}
////// creating and initializing a new tip //////////////////////
void add_newtip(int ii,int new_tip,double basex,double basey,double dir){
  double initial_tiplength;
  initial_tiplength=Min_tiplength;
  int n_pts=(int) floor(round(100000*initial_tiplength/Tip_intv)/100000);
  Tipl[new_tip]=n_pts;
  for(int j=1;j<=n_pts;j++){
    Tipx[new_tip][j]=basex+(double)(j-1)*Tip_intv*cos(dir);
    Tipy[new_tip][j]=basey+(double)(j-1)*Tip_intv*sin(dir);
  }
  Tiplength[new_tip]=calculate_curvelength(new_tip,Tipl[new_tip]);
  Tiptype[new_tip]=1;
  state[new_tip]=1;
  v[new_tip]=Vg;
  bottomneighbour[new_tip]=1;
  topneighbour[new_tip]=0;
}
///////////////////////////////////////////////////////////////////////////
///////////////////////  Branching event //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void Branching(){
  double dir,Br_dir,BrLength,bdis,prob,tiplength;
  int ii,jj,mm,ind1,ind2,ind3,base_type,old_tiplength,top1,top2;
  new_tip=n_tip;
  for (int i=1; i<=n_tip; i++){
    ////////////////////////////////////////////////////////////
    //ii=random_int(1,n_tip);///choose a random iith branch
    ii=i;
    Tiplength[ii]=calculate_curvelength(ii,Tipl[ii]);
    BrLength=BranchingLength_Threshold;
    if(Tiplength[ii]>=BrLength){///if the branch is greater than a threshold length add a new tip
      base_type=Tiptype[ii];////type of that branch
     if(ii<=N_TipInitial){
       if(Tiplength[ii]>R_Soma){
        jj=random_int((int)ceil(R_Soma/Tip_intv),Tipl[ii]-1);
        tiplength=Tiplength[ii]-R_Soma;
      }else{tiplength=0.0;jj=1;}
    }else{
      jj=random_int(3,Tipl[ii]-3);///random branching point
      tiplength=Tiplength[ii];
    }
    if(distance2D(Tipx[ii][jj],Tipy[ii][jj],0.0,0.0)>R_Soma){//branching point outside soma
      if(ran2(&idum)<1.0-exp(-tiplength*BranchingRate*Dt)){
        ///branching from  a tip, 1 new tip ,1 old tip and part of tip becomes a branch
        old_tiplength=Tipl[ii];
        //////////////////////////bottom part of tip becomes a branch //////////
        ind1=ii;ind2=new_tip+1;ind3=new_tip+2;
        Tiptype[ind1]=0;
        Tipl[ind1]=jj;
        Tiplength[ind1]=calculate_curvelength(ind1,jj);
        if(base_type==0){top1=topneighbourindex[ind1][1];top2=topneighbourindex[ind1][2];}
        topneighbour[ind1]=1;
        if(ind1>N_TipInitial){bottomneighbour[ind1]=1;}
        allocate_top_neighbour(ind1,ind2,ind3);
        //////////////////////////////upper part of tip becomes a new tip/branch
        mm=0;
        for (int j=jj; j<=old_tiplength; j++){
          mm++;
          Tipx[ind2][mm]=Tipx[ii][j];
          Tipy[ind2][mm]=Tipy[ii][j];
        }
        Tipl[ind2]=mm;
        Tiptype[ind2]=base_type;/////////new tip is either a tip or branch
        Tiplength[ind2]=calculate_curvelength(ind2,Tipl[ind2]);
        bottomneighbour[ind2]=1;
        topneighbour[ind2]=1;
        if(Tiptype[ind2]==1){//////////////tip is dynamic
          state[ind2]=state[ind1];
          v[ind2]=v[ind1];
          topneighbour[ind2]=0;
        }
        ///////////////////assign the neighbours ////////////////////////////////////
        allocate_bottom_neighbour(ind2,ind1,ind3);
        /////////////// if it is a branch reassign the neighbors for the top branch
        if(topneighbour[ind2]==1){
          allocate_top_neighbour(ind2,top1,top2);
          allocate_bottom_neighbour(top1,ind2,top2);
          allocate_bottom_neighbour(top2,ind2,top1);
        }
        //////////////////////// fresh new tip ///////////////////////////////////////
        Br_dir=wrap_angle(atan2(Tipy[ind2][2]-Tipy[ind2][1],Tipx[ind2][2]-Tipx[ind2][1]));
        //dir=wrap_angle(Br_dir+random_sign()*sample_branching_angle(Br_dir));
        dir=sample_branching_angle(Br_dir);
        add_newtip(ii,ind3,Tipx[ind2][1],Tipy[ind2][1],dir);
        //////////////////////  assign neighbours  /////////////////////////////////
        allocate_bottom_neighbour(ind3,ind1,ind2);
        ////////////increase the tip number, one branching event adds 2 more branches
        new_tip=new_tip+2;
        //////////////////////////////////////////////////////////////////////////
      }
    }
  }
}
n_tip=new_tip;
}
