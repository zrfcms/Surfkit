/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#ifdef SPaMD_viewer
#include"./../SPaMD_vis.h"
#else
#include"QB.h"
void SV_wait(char input[1024]){return;}
void SV_wait_update(double input){return;}
void SV_wait_end(){return;}
#endif
int QB_read_lmc(QB_tools *QB,const char name[50])
{
	QB_free_elename(QB);
    char gc[1024];
    FILE* input;
	int type_flag=1;
    const char id[3]="id";
    const char item[6]="ITEM:";
    const char te[5]="type";
	const char xy[3]="xy";
    const char xc[2]="x";
    const char yc[2]="y";
    const char zc[2]="z";
    const char ppc[3]="pp";
    int i,j,k;
    int extra=0;
    QB->exchannel=0;
    struct {
        int id;
        int type;
        int x;
        int y;
        int z;
    }neededdata;
	for(int i=0;i<3;i++)for(int j=0;j<3;j++)QB->mat[i][j]=0;
    neededdata.id=neededdata.type=neededdata.x=neededdata.y=neededdata.z=-1;
    input=QB_fopen(name,"r");
    int checker=0;
    do
    	if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    do
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    //set atom number
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->TotalNumber=(int)QB_checkdat(gc);
    QB->atom=(atomdata*)malloc(QB->TotalNumber*sizeof(atomdata));
	if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
	if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
	if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
	if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    //using xy yz xz pp pp pp judging box environment
    if(QB_checkchar(gc)==0)
    {
        QB->px=QB->py=QB->pz=1;
        QB->box=0;
        goto no_label;//set data when there is no label in LAMMPS output file
    }
    if(strcmp(gc,xy)==0)
    {
        QB->box=1;
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
		if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
		if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    }
    else
        QB->box=0;
    if(strcmp(gc,ppc)==0)
        QB->px=1;
    else
        QB->px=0;
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    if(strcmp(gc,ppc)==0)
        QB->py=1;
    else
        QB->py=0;
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    if(strcmp(gc,ppc)==0)
        QB->pz=1;
    else
        QB->pz=0;
    //set boundary
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}

    no_label:
    QB->startx=QB_checkdat(gc);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->boundx=QB_checkdat(gc);
    QB->boundx-=QB->startx;
    if(QB->box)
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->mat[1][0]=QB_checkdat(gc);

    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->starty=QB_checkdat(gc);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->boundy=QB_checkdat(gc);
    QB->boundy-=QB->starty;
    if(QB->box)
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->mat[2][0]=QB_checkdat(gc);

    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->startz=QB_checkdat(gc);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->boundz=QB_checkdat(gc);
    QB->boundz-=QB->startz;
    if(QB->box)
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    QB->mat[2][1]=QB_checkdat(gc);

	QB->zerox=QB->startx;
	QB->zeroy=QB->starty;
	QB->zeroz=QB->startz;
    //init lattice matrix
	if(QB->box)
	{
		QB->mat[0][0]=fabs(QB->boundx-fabs(QB->mat[1][0])-fabs(QB->mat[2][0]));
		QB->mat[1][1]=fabs(QB->boundy-fabs(QB->mat[2][1]));
		QB->mat[2][2]=QB->boundz;
		if(QB->mat[1][0]<0)
			QB->zerox-=QB->mat[1][0];
		if(QB->mat[2][0]<0)
			QB->zerox-=QB->mat[2][0];
		if(QB->mat[2][1]<0)
			QB->zeroy-=QB->mat[2][1];
	}
	else
	{
		QB->mat[0][0]=QB->boundx;
		QB->mat[1][1]=QB->boundy;
		QB->mat[2][2]=QB->boundz;
		QB->mat[0][1]=QB->mat[0][2]=QB->mat[1][0]=QB->mat[1][2]=QB->mat[2][0]=QB->mat[2][1]=0;
		QB->zerox=QB->startx;
		QB->zeroy=QB->starty;
		QB->zeroz=QB->startz;
	}
	///====================================
    ///check index data
    ///broadcast error if data doesn't start with id type x y z
    ///without those data AACSD will default to all periodic system
    ///====================================
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    //get extra data types
    j=0;
    while(1)
    {
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
        if(QB_checkchar(gc)==1)
        {
            if(strcmp(gc,id)==0)
                neededdata.id=extra;
            else if(strcmp(gc,te)==0)
				neededdata.type=extra;
            else if(strcmp(gc,xc)==0)
                neededdata.x=extra;
            else if(strcmp(gc,yc)==0)
                neededdata.y=extra;
            else if(strcmp(gc,zc)==0)
                neededdata.z=extra;
            else
            {
                printf(">>extra data %s is found\n",gc);
                QB->exchannel++;
            }
            extra++;
        }
        else
            break;
    }
    
    if(neededdata.x==-1||neededdata.y==-1||neededdata.z==-1)
    {
        printf(">>ERROR: LAMMPS dump file: %s doesn't contain enough data",name);
        exit(0);
    }
    
    QB->TypeNumber=0;
    //save data from extra channel for over write mode

    QB_fclose(input);
    input=QB_fopen(name,"r");
    do
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    do
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    do
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    do
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    while(strcmp(gc,item)!=0);
    if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
    
    QB->exch=(char**)malloc(QB->exchannel*sizeof(char*));
    for(i=0;i<QB->exchannel;i++)
        QB->exch[i]=(char*)malloc(50*sizeof(char));
    j=0;
    for(i=0;i<extra;i++)
    {
        if(QB_fscanf(input,"%s",gc)!=1){QB_fclose(input);return -1;}
        if(strcmp(gc,id)!=0&&
            strcmp(gc,te)!=0&&
            strcmp(gc,xc)!=0&&
            strcmp(gc,yc)!=0&&
            strcmp(gc,zc)!=0)
            {
                strcpy(QB->exch[j],gc); 
                j++;
            }
    }
    QB->exdat=(double**)malloc(QB->exchannel*sizeof(double*));
    for(i=0;i<QB->exchannel;i++)
        QB->exdat[i]=(double*)malloc(QB->TotalNumber*sizeof(double));
	SV_wait("Reading model please wait");
   //save needed data
    for(i=0;i<QB->TotalNumber;i++)
    {
		if((i%100000)==0)
			SV_wait_update(i*100/QB->TotalNumber);
        k=0;
        for(j=0;j<extra;j++)
        {
            if(QB_fscanf(input,"%s",gc)!=1){SV_wait_end();QB_fclose(input);return -2;}
            if(neededdata.id==j)
                QB->atom[i].id=(int)QB_checkdat(gc);
            else if(neededdata.type==j)
            {
				QB->atom[i].type=(int)QB_checkdat(gc);
				if(QB->atom[i].type<type_flag)type_flag=QB->atom[i].type;
			}
            else if(neededdata.x==j)
                QB->atom[i].x=QB_checkdat(gc);
            else if(neededdata.y==j)
                QB->atom[i].y=QB_checkdat(gc);
            else if(neededdata.z==j)
                QB->atom[i].z=QB_checkdat(gc);
            else if(QB->exchannel!=0)
            {
                QB->exdat[k][i]=QB_checkdat(gc);
                k++;
            }  
            if(neededdata.id==-1)
                QB->atom[i].id=i+1;
            if(neededdata.type==-1)
                QB->atom[i].type=1;
        }
        if(QB->atom[i].type>QB->TypeNumber)
            QB->TypeNumber=QB->atom[i].type;
    }
	for(i=0;i<QB->TotalNumber;i++)
    {
		QB->atom[i].type=QB->atom[i].type+1-type_flag;
	}
	QB->TypeNumber=QB->TypeNumber+1-type_flag;
    QB_printf(">>LAMMPS dump file %s is found\n",name);
    QB_printf(">>total %d atoms\n",QB->TotalNumber);
    QB_printf(">>total %d atom types\n",QB->TypeNumber);
    QB_fclose(input);
	SV_wait_end();
    if(QB_checkdat_err){QB_checkdat_err=0;return -3;}//char data in numbers
    return 1;
}
//==========================================================================
void QB_free_exchannel(QB_tools* QB)
{
    int i;
    if(QB->exchannel)
    {
        for(i=0;i<QB->exchannel;i++)
            free(QB->exdat[i]);
        for(i=0;i<QB->exchannel;i++)
            free(QB->exch[i]);
        free(QB->exdat);
        free(QB->exch);
        QB->exchannel=0;
    }
}
