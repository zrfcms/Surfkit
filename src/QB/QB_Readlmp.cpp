/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#ifdef SPaMD_viewer
#include"./../SPaMD_vis.h"
#else
#include"QB.h"
extern void SV_wait(char input[1024]);
extern void SV_wait_update(double input);
extern void SV_wait_end();
#endif
int QB_read_lmp(QB_tools *QB,const char name[50])
{    
	QB_free_elename(QB);
    FILE    *fi2;
    double  x,y,z;
    int     n,now;
    char    ch[256], cc[6]="Atoms";
    QB->px=QB->py=QB->pz=1;
	for(int i=0;i<3;i++)for(int j=0;j<3;j++)QB->mat[i][j]=0;
    fi2=QB_fopen(name,"r");
    QB_fscanf(fi2,"%*[^\n]%*c");
    if(QB_fscanf(fi2,"%d ",&QB->TotalNumber)!=1){QB_fclose(fi2);return -1;}
    QB->atom=(atomdata *)malloc(QB->TotalNumber*sizeof(atomdata));
    QB_fscanf(fi2,"%*[^\n]%*c");
    if(QB_fscanf(fi2,"%d ",&QB->TypeNumber)!=1){QB_fclose(fi2);return -1;}
    QB_fscanf(fi2,"%*[^\n]%*c");
    if(QB_fscanf(fi2,"%lf%lf ",&QB->startx,&QB->boundx)!=2){QB_fclose(fi2);return -1;}
    QB->boundx=QB->boundx-QB->startx;
    QB_fscanf(fi2,"%*[^\n]%*c");
    if(QB_fscanf(fi2,"%lf%lf ",&QB->starty,&QB->boundy)!=2){QB_fclose(fi2);return -1;}
    QB->boundy=QB->boundy-QB->starty;
    QB_fscanf(fi2,"%*[^\n]%*c");
    if(QB_fscanf(fi2,"%lf%lf ",&QB->startz,&QB->boundz)!=2){QB_fclose(fi2);return -1;}
    QB->boundz=QB->boundz-QB->startz;
    fscanf(fi2,"%*[^\n]%*c");
    QB->box=QB_fscanf(fi2,"%lf%lf%lf ",&(QB->mat[1][0]),&(QB->mat[2][0]),&(QB->mat[2][1]))/3;

    while(1)
    {
        if(QB_fscanf(fi2,"%s",ch)!=1){QB_fclose(fi2);return -1;}
        if(strcmp(ch,cc)==0)
            break;
    }
	SV_wait("Reading model please wait");
    QB_fscanf(fi2,"%*[^\n]%*c");
    for(now=0;now<QB->TotalNumber;now++)
    {
		if((now%100000)==0)
			SV_wait_update(now*100/QB->TotalNumber);
        if(QB_fscanf(fi2,"%d%u%lf%lf%lf ",&QB->atom[now].id,&n,&x,&y,&z)!=5){SV_wait_end();QB_fclose(fi2);return -2;}
        QB->atom[now].x=x;
        QB->atom[now].y=y;
        QB->atom[now].z=z;
        QB->atom[now].type=n;
    }
	QB->mat[0][0]=QB->boundx;
    QB->mat[1][1]=QB->boundy;
    QB->mat[2][2]=QB->boundz;
	if(QB->box)
	{
		QB->zerox=QB->startx;
		QB->zeroy=QB->starty;
		QB->zeroz=QB->startz;
		if(QB->mat[1][0]<0)QB->startx+=QB->mat[1][0];
		if(QB->mat[2][0]<0)QB->startx+=QB->mat[2][0];
		if(QB->mat[2][1]<0)QB->starty+=QB->mat[2][1];
		QB->boundx+=fabs(QB->mat[1][0])+fabs(QB->mat[2][0]);
		QB->boundy+=fabs(QB->mat[2][1]);
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
    QB_printf(">>LAMMPS dump file %s is found\n",name);
    QB_printf(">>total %d atoms\n",QB->TotalNumber);
    QB_printf(">>total %d atom types\n",QB->TypeNumber);
    QB_fclose(fi2);
	SV_wait_end();
    return 1;
}
