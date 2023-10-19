void removeSpaces(char *str){
    int count = 0;
    for (int i = 0; str[i]; i++)
        if (str[i] != ' ')
            str[count++] = str[i]; /// incremented
    str[count] = '\0';
}

int ReadInputData(int argc,char *arg0, char *fname) {
if(argc<2){
  fprintf(stderr,"Usage %s Parameters.in\n",arg0);
  return -1;
}
FILE *fp;
fp=fopen(fname,"r");;
if(fp==NULL){
  printf("Unable to open file %s\n",fname);
  return -1;
}
fclose(fp);
/////////////////////// Read data /////////////////////////////////
char *str[]={"N_Sample","Dump_Conf","Time_Steps","Dump_Data",
"Tip_Persis","BranchingAngleMean","BranchingAngleStd","BranchingRate","Vg","Vs","Kgp","Kgs","Kpg","Kps","Ksg","Ksp","SizeY","SizeX","pixelsize",
"Boundary_type","Print_Conf"};
int n_parameters=21;
    int vcount;
	int line_num ;
	int find_result;
    char* saveptr = NULL;
    char temp[512];
    char *varval;
    const char delimiters[] = ",;:!-=";
    char *token; 
for(int n=0;n<n_parameters;n++){
    fp=fopen(fname,"r");;
    line_num=1;
    find_result=0;
	while(fgets(temp, 512, fp) != NULL) { 
		if((strstr(temp, str[n])) != NULL) {
            token = strtok_r(temp, delimiters, &saveptr);
            vcount=0;
            while (token != NULL) {
                vcount++;
                if(vcount==2){
                    removeSpaces(token);
                    if(n==0){N_Sample=atoi(token);}
                    if(n==1){Dump_Conf=atoi(token);}
                    if(n==2){Time_Steps=atoi(token);}
                    if(n==3){Dump_Data=atoi(token);}
                    if(n==4){Tip_Persis=atof(token);}
                    if(n==5){BranchingAngleMean=atof(token)*pi/180.0;;}
                    if(n==6){BranchingAngleStd=atof(token)*pi/180;}
                    if(n==7){BranchingRate=atof(token);}
                    if(n==8){Vg=atof(token);}
                    if(n==9){Vs=fabs(atof(token));}
                    if(n==10){Kgp=atof(token);}
                    if(n==11){Kgs=atof(token);}
                    if(n==12){Kpg=atof(token);}
                    if(n==13){Kps=atof(token);}
                    if(n==14){Ksg=atof(token);}
                    if(n==15){Ksp=atof(token);}
                    if(n==16){SizeY=atof(token);}
                    if(n==17){SizeX=atof(token);}
                    if(n==18){pixelsize=atof(token);}
                    if(n==19){strcpy(Boundary_type,token);}
                    if(n==20){strcpy(Print_Conf,token);}
                    }
                token = strtok_r(NULL,delimiters,&saveptr);
            }
		find_result++;
		}
		line_num++;
	}
	if(find_result == 0) {printf("\nCouldn't find parameter %s, Please enter it in this format\n%s=value;\n",str[n],str[n]);return -1;}
    fclose(fp);
}
if(fp) {fclose(fp);}//Close the file if still open
return(0);
}
///////////////////////////////////////////
