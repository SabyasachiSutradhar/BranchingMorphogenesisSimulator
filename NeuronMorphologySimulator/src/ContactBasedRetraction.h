/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022, as a part of
Neuonal morphogenesis main code.
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
This header excutes the non ovelaping interaction among the branches and the boundary.
Copyright @ Sabyasachi Sutradhar
*/

double Return_Radius(int i){
  /*
  double rr;
  if(i<=N_TipInitial){
    rr=0.5;
  }else{
  rr=0.1+Tipage[i]*0.2/(100.0*60.0);
  if(rr>0.3){rr=0.3;}
}
*/
  return R_Branch;
}

double pbc(double x,double Lx){
  if(strcmp(Boundary_type,"pbc")==0){
    if(x>Lx){x=x-2.0*Lx;}
    if(x<-Lx){x=2.0*Lx+x;}
}
  return x;
}

double distance_pbc(double x1,double y1,double x2,double y2){
double dis,dx,dy;
//dx=(pbc(x1,SizeAP)-pbc(x2,SizeAP));
//dy=(pbc(y1,SizeLR)-pbc(y2,SizeLR));
/////pbc
dx=x1-x2;
dy=y1-y2;

if(strcmp(Boundary_type,"pbc")==0){
  if (fabs(dx)>SizeX/2.0){dx=2.0*SizeX-fabs(dx);}
  if (fabs(dy)>SizeY/2.0){dy=2.0*SizeY-fabs(dy);}
}
dis=sqrt(dx*dx+dy*dy);
return dis;
}

int pbc_on(double x1,double y1,double x2,double y2){
  int ind=0;
  double dx,dy;
  dx=(x1-x2);
  dy=(y1-y2);
    if (fabs(dx)>SizeX || fabs(dy)>SizeY ){ind=1;}
  return ind;
}
////////////////////////////////////////////////////////////////////////////////////////////
  double Dis_Line2Point(double x1,double y1,double x2,double y2,double x0,double y0){
  return fabs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
  }
////////////////////////////////////////////////////////////////////////////////////////////
///// determine whether tip hits other tips//////////////////////////////////////////////////
  int branch_collision(int i,double xt,double yt){
    int ind=0,j1,j2;
    int ne;
    double cmx,cmy,ri,rj,xl,xh,yl,yh,L_box;
    ri=Return_Radius(i);
///////////////////////////////////////// self interaction ///////////////////////////////////////
        for(int k=1;k<Tipl[i]-1;k++){
          if(k>=1 && k<=Tipl[i]-10){
          if(distance2D(Tipx[i][k],Tipy[i][k],xt,yt)<=Interaction_length){ind=1;break;}
          if(lineseg_intersection(Tipx[i][k],Tipy[i][k],Tipx[i][k+1],Tipy[i][k+1],Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],xt,yt)==1){ind=1;break;}
      }
    }
///////////////// interaction with the neighbours//////////////////////////////////
      j1=bottomneighbourindex[i][1]; j2=bottomneighbourindex[i][2];

      for(int k=1;k<Tipl[j1];k++){
        if(distance2D(Tipx[j1][k],Tipy[j1][k],xt,yt)<=Interaction_length){ind=1;break;}
        if(lineseg_intersection(Tipx[j1][k],Tipy[j1][k],Tipx[j1][k+1],Tipy[j1][k+1],Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],xt,yt)==1){ind=1;break;}
      }
      for(int k=1;k<Tipl[j2];k++){
    if(distance2D(Tipx[j2][k],Tipy[j2][k],xt,yt)<=Interaction_length){ind=1;break;}
        if(lineseg_intersection(Tipx[j2][k],Tipy[j2][k],Tipx[j2][k+1],Tipy[j2][k+1],Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],xt,yt)==1){ind=1;break;}
      }
//////////////////// all other branch within certain radius to minimize the computation time
//  if((int)(time/Dt)%2==0){
    for(int j=1;j<=n_tip;j++){
      if(j!=i && j!=j1 && j!=j2){
        cmx=(Tipx[j][1]+Tipx[j][Tipl[j]])/2.0;
        cmy=(Tipy[j][1]+Tipy[j][Tipl[j]])/2.0;

        if(distance_pbc(Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],cmx,cmy)<=25.0 || distance_pbc(Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],Tipx[j][1],Tipy[j][1])<=25.0 || distance_pbc(Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],Tipx[j][Tipl[j]],Tipy[j][Tipl[j]])<=25.0){
          rj=Return_Radius(j);
          for(int k=1;k<Tipl[j];k++){
            if(distance_pbc(Tipx[j][k],Tipy[j][k],xt,yt)<=ri+rj+Interaction_length){ind=1;break;}
            if(lineseg_intersection(Tipx[j][k],Tipy[j][k],Tipx[j][k+1],Tipy[j][k+1],Tipx[i][Tipl[i]],Tipy[i][Tipl[i]],xt,yt)==1){ind=1;break;}
            //if(lineseg_intersection(pbc(Tipx[j][k],SizeAP),pbc(Tipy[j][k],SizeLR),pbc(Tipx[j][k+1],SizeAP),pbc(Tipy[j][k+1],SizeLR),pbc(Tipx[i][Tipl[i]],SizeAP),pbc(Tipy[i][Tipl[i]],SizeLR),pbc(xt,SizeAP),pbc(yt,SizeLR))==1){ind=1;break;}
            }
          }
        }
      }
     //}
  return ind;
  }


/////// determine if the tip hits the soma///////////////////////////////////////////////////
  int soma_collision(int i,double xt,double yt){
    int ind=0;
    if(Tiptype[i]==1){
      if(distance2D(0.0,0.0,xt,yt)<=R_Soma){
        ind=1;
    }
  }
    return ind;
}
