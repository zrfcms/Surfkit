#include"./../QB/QB.h"
#include"QSPG.h"

void QB_c2f(QB_tools *QB,double input[3],double output[3])
{
	double temp[3];
	double A;
	double mat[3][3];
	if(QB->box)
	{
		A=(QB->mat[0][0]*QB->mat[1][1]*QB->mat[2][2]+QB->mat[0][1]*QB->mat[1][2]*QB->mat[2][0]+QB->mat[0][2]*QB->mat[1][0]*QB->mat[2][1]
		  -QB->mat[0][0]*QB->mat[1][2]*QB->mat[2][1]-QB->mat[0][1]*QB->mat[1][0]*QB->mat[2][2]-QB->mat[0][2]*QB->mat[1][1]*QB->mat[2][0]);
		if(A==0){output[0]=output[1]=output[2]=0;return;}
		A=1/A;
		mat[0][0]=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
		mat[0][1]=-QB->mat[1][0]*QB->mat[2][2]+QB->mat[1][2]*QB->mat[2][0];
		mat[0][2]=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
		mat[1][0]=-QB->mat[0][1]*QB->mat[2][2]+QB->mat[0][2]*QB->mat[2][1];
		mat[1][1]=QB->mat[0][0]*QB->mat[2][2]-QB->mat[0][2]*QB->mat[2][0];
		mat[1][2]=-QB->mat[0][0]*QB->mat[2][1]+QB->mat[0][1]*QB->mat[2][0];
		mat[2][0]=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
		mat[2][1]=-QB->mat[0][0]*QB->mat[1][2]+QB->mat[0][2]*QB->mat[1][0];
		mat[2][2]=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
	}
	else
	{
		temp[0]=1.0/QB->boundx;
		temp[1]=1.0/QB->boundy;
		temp[2]=1.0/QB->boundz;
	}
	if(QB->box)
	{
		output[0]=A*(input[0]*mat[0][0]+input[1]*mat[0][1]+input[2]*mat[0][2]);
		output[1]=A*(input[0]*mat[1][0]+input[1]*mat[1][1]+input[2]*mat[1][2]);
		output[2]=A*(input[0]*mat[2][0]+input[1]*mat[2][1]+input[2]*mat[2][2]);
	}
	else
	{
		output[0]=temp[0]*input[0];
		output[1]=temp[1]*input[1];
		output[2]=temp[2]*input[2];
	}
}

void QB_f2c(QB_tools *QB,double input[3],double output[3])
{
	if(QB->box)
	{
		output[0]=input[0]*QB->mat[0][0]+input[1]*QB->mat[1][0]+input[2]*QB->mat[2][0]+QB->zerox;
		output[1]=input[0]*QB->mat[0][1]+input[1]*QB->mat[1][1]+input[2]*QB->mat[2][1]+QB->zeroy;
		output[2]=input[0]*QB->mat[0][2]+input[1]*QB->mat[1][2]+input[2]*QB->mat[2][2]+QB->zeroz;
	}
	else
	{
		output[0]=input[0]*QB->boundx+QB->startx;
		output[1]=input[1]*QB->boundy+QB->starty;
		output[2]=input[2]*QB->boundz+QB->startz;
	}
}

void QB2SPG(QB_tools *QB,SPG_tools* SPG)
{
	int i,j,k;
	double temp[3];
	double A;
	double mat[3][3];
	
	if(SPG->num!=0)
	{
		free(SPG->pos);
		free(SPG->type);
	}
	SPG->num=0;
	
	if(SPG->ele_n!=0)
	{
		free(SPG->elename);
	}
	SPG->ele_n=0;
	
	for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		SPG->mat[j][i]=QB->mat[i][j];
		
	if(QB->TotalNumber)
	{
		SPG->pos=(double(*)[3])malloc(sizeof(double)*3*QB->TotalNumber);
		SPG->type=(int*)malloc(sizeof(int)*QB->TotalNumber);
	}
	if(QB->TypeNumber)
		SPG->elename=(char(*)[32])malloc(sizeof(char)*32*QB->TypeNumber);
	SPG->ele_n=QB->TypeNumber;
	for(i=0;i<QB->TypeNumber;i++)
	{
		QB_load_elename(QB,i,SPG->elename[i]);
	}
	
	if(QB->box)
	{
		A=1.00/(QB->mat[0][0]*QB->mat[1][1]*QB->mat[2][2]+QB->mat[0][1]*QB->mat[1][2]*QB->mat[2][0]+QB->mat[0][2]*QB->mat[1][0]*QB->mat[2][1]
		       -QB->mat[0][0]*QB->mat[1][2]*QB->mat[2][1]-QB->mat[0][1]*QB->mat[1][0]*QB->mat[2][2]-QB->mat[0][2]*QB->mat[1][1]*QB->mat[2][0]);
		mat[0][0]=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
		mat[0][1]=-QB->mat[1][0]*QB->mat[2][2]+QB->mat[1][2]*QB->mat[2][0];
		mat[0][2]=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
		mat[1][0]=-QB->mat[0][1]*QB->mat[2][2]+QB->mat[0][2]*QB->mat[2][1];
		mat[1][1]=QB->mat[0][0]*QB->mat[2][2]-QB->mat[0][2]*QB->mat[2][0];
		mat[1][2]=-QB->mat[0][0]*QB->mat[2][1]+QB->mat[0][1]*QB->mat[2][0];
		mat[2][0]=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
		mat[2][1]=-QB->mat[0][0]*QB->mat[1][2]+QB->mat[0][2]*QB->mat[1][0];
		mat[2][2]=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
	}
	else
	{
		temp[0]=1.0/QB->boundx;
		temp[1]=1.0/QB->boundy;
		temp[2]=1.0/QB->boundz;
	}
	k=0;
    for(i=0;i<QB->TotalNumber;i++)
    {
        if(QB->atom[i].id>0)
		{
			SPG->type[k]=QB->atom[i].type;
			if(QB->box)
			{
				SPG->pos[k][0]=A*((QB->atom[i].x-QB->zerox)*mat[0][0]+(QB->atom[i].y-QB->zeroy)*mat[0][1]+(QB->atom[i].z-QB->zeroz)*mat[0][2]);
				SPG->pos[k][1]=A*((QB->atom[i].x-QB->zerox)*mat[1][0]+(QB->atom[i].y-QB->zeroy)*mat[1][1]+(QB->atom[i].z-QB->zeroz)*mat[1][2]);
				SPG->pos[k][2]=A*((QB->atom[i].x-QB->zerox)*mat[2][0]+(QB->atom[i].y-QB->zeroy)*mat[2][1]+(QB->atom[i].z-QB->zeroz)*mat[2][2]);
			}
			else
			{
				SPG->pos[k][0]=temp[0]*(QB->atom[i].x-QB->startx);
				SPG->pos[k][1]=temp[1]*(QB->atom[i].y-QB->starty);
				SPG->pos[k][2]=temp[2]*(QB->atom[i].z-QB->startz);
				
			}
			k++;
		}
    }
	SPG->num=k;
}

void SPG2QB(SPG_tools* SPG,QB_tools *QB)
{
	int i,j,k;
	double temp[3];
	double A;
	double mat[3][3];
	int typelist_n=0;
	if(QB->TotalNumber!=0)
	{
		QB_free_atom(QB);
	}
	QB->TotalNumber=SPG->num;
	for(i=0;i<3;i++)for(j=0;j<3;j++)QB->mat[i][j]=SPG->mat[j][i];
	QB->px=QB->py=QB->pz=1;
	QB->box=1;
	QB->boundx=fabs(QB->mat[0][0])+fabs(QB->mat[1][0])+fabs(QB->mat[2][0]);
    QB->boundy=fabs(QB->mat[0][1])+fabs(QB->mat[1][1])+fabs(QB->mat[2][1]);
	QB->boundz=fabs(QB->mat[0][2])+fabs(QB->mat[1][2])+fabs(QB->mat[2][2]);
	QB->startx=QB->zerox=QB->starty=QB->zeroy=QB->startz=QB->zeroz=0;
	if(QB->mat[0][0]<0)QB->startx+=QB->mat[0][0];
	if(QB->mat[1][0]<0)QB->startx+=QB->mat[1][0];
	if(QB->mat[2][0]<0)QB->startx+=QB->mat[2][0];
	if(QB->mat[0][1]<0)QB->starty+=QB->mat[0][1];
	if(QB->mat[1][1]<0)QB->starty+=QB->mat[1][1];
	if(QB->mat[2][1]<0)QB->starty+=QB->mat[2][1];
	if(QB->mat[0][2]<0)QB->startz+=QB->mat[0][2];
	if(QB->mat[1][2]<0)QB->startz+=QB->mat[1][2];
	if(QB->mat[2][2]<0)QB->startz+=QB->mat[2][2];
	if(SPG->num)
		QB->atom=(atomdata*)malloc(SPG->num*sizeof(atomdata));

    for(i=0;i<QB->TotalNumber;i++)
    {
		QB->atom[i].id=i+1;
		QB->atom[i].type=SPG->type[i];
		if(QB->atom[i].type>typelist_n)typelist_n=QB->atom[i].type;
		QB->atom[i].x=SPG->pos[i][0]*QB->mat[0][0]+SPG->pos[i][1]*QB->mat[1][0]+SPG->pos[i][2]*QB->mat[2][0];
		QB->atom[i].y=SPG->pos[i][0]*QB->mat[0][1]+SPG->pos[i][1]*QB->mat[1][1]+SPG->pos[i][2]*QB->mat[2][1];
		QB->atom[i].z=SPG->pos[i][0]*QB->mat[0][2]+SPG->pos[i][1]*QB->mat[1][2]+SPG->pos[i][2]*QB->mat[2][2];
    }
	QB->TypeNumber=typelist_n;
	
	for(i=0;i<SPG->ele_n;i++)
	{
		QB_save_elename(QB,i,SPG->elename[i]);
	}
}

void SPG_free(SPG_tools* SPG)
{
	if(SPG->num!=0)
	{
		free(SPG->pos);
		free(SPG->type);
	}
	SPG->num=0;
}