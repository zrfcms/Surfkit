/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
#ifdef SPaMD_model
#include"./../../SPaMD_basic/SPaMD_globel.h"
//0.031685 M1
#endif
void QB_init(QB_tools *QB)
{
	QB->num_exint=QB->num_exdouble=QB->TotalNumber=QB->exchannel=QB->TypeNumber=0;
	QB->MCN=16;
	QB->box=0;
	QB->mat[0][0]=QB->mat[0][1]=QB->mat[0][2]=QB->mat[1][0]=QB->mat[1][1]=QB->mat[1][2]=QB->mat[2][0]=QB->mat[2][1]=QB->mat[2][2]=0;
	QB->network_flag=QB->list_flag=QB->ele.n=0;
	QB->px=QB->py=QB->pz=1;
}

double QB_get_data(QB_tools *QB,int i,const char *name)
{
	int j;
	if(strcmp(name,"id")==0)
		return (double)(QB->atom[i].id);
	if(strcmp(name,"type")==0)
		return (double)(QB->atom[i].type);
	if(strcmp(name,"x")==0)
		return QB->atom[i].x;
	if(strcmp(name,"y")==0)
		return QB->atom[i].y;
	if(strcmp(name,"z")==0)
		return QB->atom[i].z;
	for(j=0;j<QB->exchannel;j++)
		if(strcmp(name,QB->exch[j])==0)
		{
			return QB->exdat[j][i];
		}
	for(j=0;j<QB->num_exint;j++)
		if(strcmp(name,QB->name_exint[j])==0)
		{
			return (double)(QB->par_exint[j][i]);
		}
	for(j=0;j<QB->num_exdouble;j++)
		if(strcmp(name,QB->name_exdouble[j])==0)
		{
			return QB->par_exdouble[j][i];
		}
	return 0;
}

void QB_save_data(QB_tools *QB,int i,const char *name,double x)
{
	int j;
	for(j=0;j<QB->exchannel;j++)
		if(strcmp(name,QB->exch[j])==0)
		{
			QB->exdat[j][i]=x;
			return;
		}
	for(j=0;j<QB->num_exint;j++)
		if(strcmp(name,QB->name_exint[j])==0)
		{
			QB->par_exint[j][i]=(int)x;
			return;
		}
	for(j=0;j<QB->num_exdouble;j++)
		if(strcmp(name,QB->name_exdouble[j])==0)
		{
			QB->par_exdouble[j][i]=x;
			return;
		}
}
//for massive data reading
int QB_slot_sum(QB_tools *QB)
{
	return 5+QB->exchannel+QB->num_exint+QB->num_exdouble;
}

int QB_slot_get(QB_tools *QB,const char *name)
{
	if(strcmp(name,"id")==0)
		return 0;
	if(strcmp(name,"type")==0)
		return 1;
	if(strcmp(name,"x")==0)
		return 2;
	if(strcmp(name,"y")==0)
		return 3;
	if(strcmp(name,"z")==0)
		return 4;
	int j;
	for(j=0;j<QB->exchannel;j++)
		if(strcmp(name,QB->exch[j])==0)
		{
			return 5+j;
		}
	for(j=0;j<QB->num_exint;j++)
		if(strcmp(name,QB->name_exint[j])==0)
		{
			return 5+QB->exchannel+j;
		}
	for(j=0;j<QB->num_exdouble;j++)
		if(strcmp(name,QB->name_exdouble[j])==0)
		{
			return 5+QB->exchannel+QB->num_exint+j;
		}
	return -1;
}

int QB_slot_name(QB_tools *QB,int i,char *name)
{
	if(i==0){sprintf(name,"id");return 1;}
	else if(i==1){sprintf(name,"type");return 1;}
	else if(i==2){sprintf(name,"x");return 1;}
	else if(i==3){sprintf(name,"y");return 1;}
	else if(i==4){sprintf(name,"z");return 1;}
	else if(i<5+QB->exchannel){sprintf(name,"%s",QB->exch[i-5]);return 1;}
	else if(i<5+QB->exchannel+QB->num_exint){sprintf(name,"%s",QB->name_exint[i-5-QB->exchannel]);return 1;}
	else if(i<5+QB->exchannel+QB->num_exint+QB->num_exdouble){sprintf(name,"%s",QB->name_exdouble[i-5-QB->exchannel-QB->num_exint]);return 1;}
	return 0;
}

double QB_slot_data(QB_tools *QB,int j,int i)
{
	if(j==0)return QB->atom[i].id;
	else if(j==1)return QB->atom[i].type;
	else if(j==2)return QB->atom[i].x;
	else if(j==3)return QB->atom[i].y;
	else if(j==4)return QB->atom[i].z;
	else if(j<5+QB->exchannel)return QB->exdat[j-5][i];
	else if(j<5+QB->exchannel+QB->num_exint)return (double)(QB->par_exint[j-5-QB->exchannel][i]);
	else if(j<5+QB->exchannel+QB->num_exint+QB->num_exdouble)return QB->par_exdouble[j-5-QB->exchannel-QB->num_exint][i];
	return 0;
}

void QB_slot_save(QB_tools *QB,int j,int i,double x)
{
	if(j==0)QB->atom[i].id=(int)x;
	else if(j==1)QB->atom[i].type=(int)x;
	else if(j==2)QB->atom[i].x=x;
	else if(j==3)QB->atom[i].y=x;
	else if(j==4)QB->atom[i].z=x;
	else if(j<5+QB->exchannel)QB->exdat[j-5][i]=x;
	else if(j<5+QB->exchannel+QB->num_exint)QB->par_exint[j-5-QB->exchannel][i]=(int)x;
	else if(j<5+QB->exchannel+QB->num_exint+QB->num_exdouble)QB->par_exdouble[j-5-QB->exchannel-QB->num_exint][i]=x;
}

void QB_add_exint(QB_tools *QB,const char *name)
{
	if(QB->TotalNumber==0)
	{
		printf("Please get at least one number before add extra data\n");
		return;
	}
	int i;
    if(QB_slot_get(QB,name)!=-1)
        return;
	if(QB->num_exint==0)
	{
		QB->name_exint=(char**)malloc(sizeof(char*));
		QB->par_exint=(int**)malloc(sizeof(int*));
	}
	else 
	{
		QB->name_exint=(char**)realloc(QB->name_exint,(QB->num_exint+1)*sizeof(char*));
		QB->par_exint=(int**)realloc(QB->par_exint,(QB->num_exint+1)*sizeof(int*));
	}
	QB->name_exint[QB->num_exint]=(char*)malloc(50*sizeof(char));
	QB->par_exint[QB->num_exint]=(int*)malloc(QB->TotalNumber*sizeof(int));
	strcpy(QB->name_exint[QB->num_exint],name);
	QB->num_exint++;
}

void QB_add_exdouble(QB_tools *QB,const char *name)
{
	if(QB->TotalNumber==0)
	{
		QB_printf("Please get at least one atom before add extra data\n");
		return;
	}
    int i;
    if(QB_slot_get(QB,name)!=-1)
        return;
	if(QB->num_exdouble==0)
	{
		QB->name_exdouble=(char**)malloc(sizeof(char*));
		QB->par_exdouble=(double**)malloc(sizeof(double*));
	}
	else 
	{
		QB->name_exdouble=(char**)realloc(QB->name_exdouble,(QB->num_exdouble+1)*sizeof(char*));
		QB->par_exdouble=(double**)realloc(QB->par_exdouble,(QB->num_exdouble+1)*sizeof(double*));
	}
	QB->name_exdouble[QB->num_exdouble]=(char*)malloc(50*sizeof(char));
	QB->par_exdouble[QB->num_exdouble]=(double*)malloc(QB->TotalNumber*sizeof(double));
	strcpy(QB->name_exdouble[QB->num_exdouble],name);
	QB->num_exdouble++;
}

void QB_free_extra(QB_tools *QB)
{
	int i;
	for(i=0;i<QB->num_exint;i++)
	{
		free(QB->name_exint[i]);
		free(QB->par_exint[i]);
	}
	for(i=0;i<QB->num_exdouble;i++)
	{
		free(QB->name_exdouble[i]);
		free(QB->par_exdouble[i]);
	}	
	if(QB->num_exint!=0)
	{
		free(QB->name_exint);
		free(QB->par_exint);
		QB->num_exint=0;
	}
	if(QB->num_exdouble!=0)
	{
		free(QB->name_exdouble);
		free(QB->par_exdouble);
		QB->num_exdouble=0;
	}
}

void QB_free_atom(QB_tools *QB)
{
	if(QB->TotalNumber>0)
		free(QB->atom);
	QB->TotalNumber=0;
	QB_free_extra(QB);
	QB_free_exchannel(QB);
	QB_free_elename(QB);
	QB->num_exint=QB->num_exdouble=QB->TotalNumber=QB->exchannel=QB->TypeNumber=0;
	QB->network_flag=QB->list_flag=QB->ele.n=0;
	QB->px=QB->py=QB->pz=1;
}

void QB_free_all(QB_tools *QB)
{
	QB_free_atom(QB);
	QB_init(QB);
}

void QB_delete_atom(QB_tools *QB,int i)
{
	if(QB->atom[i].id>0)
		QB->atom[i].id=-QB->atom[i].id;
}

void QB_recover_atom(QB_tools *QB,int i)
{
	if(QB->atom[i].id<0)
		QB->atom[i].id=-QB->atom[i].id;
}

void QB_update_extra(QB_tools *QB)
{
	int i;
	for(i=0;i<QB->exchannel;i++)
    	        QB->exdat[i]=(double*)realloc(QB->exdat[i],QB->TotalNumber*sizeof(double));
    	for(i=0;i<QB->num_exint;i++)
    	      	QB->par_exint[i]=(int*)realloc(QB->par_exint[i],QB->TotalNumber*sizeof(int));
    	for(i=0;i<QB->num_exdouble;i++)
    	      	QB->par_exdouble[i]=(double*)realloc(QB->par_exdouble[i],QB->TotalNumber*sizeof(double));
}
#define DELTA_ADD_ATOM 1
void QB_create_atom(QB_tools *QB,int type,double x,double y,double z)
{
	int i;
	if(QB->TotalNumber==0)
	{
		QB->atom=(atomdata*)malloc(DELTA_ADD_ATOM*sizeof(atomdata));
		QB->atom[0].id=1;
		QB->atom[0].type=type;
		QB->atom[0].x=x;
		QB->atom[0].y=y;
		QB->atom[0].z=z;
		QB->TotalNumber++;
	}
	else
	{
		if(QB->TotalNumber%DELTA_ADD_ATOM==0)
		QB->atom=(atomdata*)realloc(QB->atom,(QB->TotalNumber+DELTA_ADD_ATOM)*sizeof(atomdata));
		QB->atom[QB->TotalNumber].id=QB->TotalNumber+1;
		QB->atom[QB->TotalNumber].type=type;
		QB->atom[QB->TotalNumber].x=x;
		QB->atom[QB->TotalNumber].y=y;
		QB->atom[QB->TotalNumber].z=z;
		for(i=0;i<QB->exchannel;i++)
    	        	QB->exdat[i]=(double*)realloc(QB->exdat[i],(QB->TotalNumber+1)*sizeof(double));
    	        for(i=0;i<QB->num_exint;i++)
    	        	QB->par_exint[i]=(int*)realloc(QB->par_exint[i],(QB->TotalNumber+1)*sizeof(int));
    	        for(i=0;i<QB->num_exdouble;i++)
    	        	QB->par_exdouble[i]=(double*)realloc(QB->par_exdouble[i],(QB->TotalNumber+1)*sizeof(double));
    	        QB->TotalNumber++;
	}	
	if(type>QB->TypeNumber)
            QB->TypeNumber=type;
}

void QB_swap(QB_tools *QB,int i,int j)
{
	int k;
	int temp_int;
	double temp_double;
	
	temp_int=QB->atom[i].id;
	QB->atom[i].id=QB->atom[j].id;
	QB->atom[j].id=temp_int;
	temp_int=QB->atom[i].type;
	QB->atom[i].type=QB->atom[j].type;
	QB->atom[j].type=temp_int;
	temp_double=QB->atom[i].x;
	QB->atom[i].x=QB->atom[j].x;
	QB->atom[j].x=temp_double;
	temp_double=QB->atom[i].y;
	QB->atom[i].y=QB->atom[j].y;
	QB->atom[j].y=temp_double;
	temp_double=QB->atom[i].z;
	QB->atom[i].z=QB->atom[j].z;
	QB->atom[j].z=temp_double;
	for(k=0;k<QB->exchannel;k++)
    	{
    		temp_double=QB->exdat[k][i];
    		QB->exdat[k][i]=QB->exdat[k][j];
    		QB->exdat[k][j]=temp_double;
    	}
    	for(k=0;k<QB->num_exint;k++)
    	{
    		temp_int=QB->par_exint[k][i];
    		QB->par_exint[k][i]=QB->par_exint[k][j];
    		QB->par_exint[k][j]=temp_int;
    	}
    	for(k=0;k<QB->num_exdouble;k++)
    	{
    		temp_double=QB->par_exdouble[k][i];
    		QB->par_exdouble[k][i]=QB->par_exdouble[k][j];
    		QB->par_exdouble[k][j]=temp_double;
    	}
}


void QB_affine_delete(QB_tools *QB)
{
	int i,j;
	int initial_id=QB->TotalNumber;
	for(i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].id<0)
			for(j=QB->TotalNumber-1;j>=0;j--)
			{
				if(i==j)
				{	
					QB->TotalNumber--;
					break;
				}
				if(QB->atom[j].id<0)
				{
					QB->TotalNumber--;
					continue;
				}
				else
				{
					QB_swap(QB,i,j);
					QB->TotalNumber--;
					break;
				}
			}
	}
	if(initial_id!=QB->TotalNumber)
	{
		for(i=0;i<QB->TotalNumber;i++)
		{
			QB->atom[i].id=i+1;
		}
		QB->atom=(atomdata*)realloc(QB->atom,QB->TotalNumber*sizeof(atomdata));
	}
	
}

void QB_dump_lmp(QB_tools *QB,const char *name)
{
	int i,all=0,k=0;
	for(i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].id>0)
			all++;
	}
	FILE* outfile;
	outfile=QB_fopen(name,"w");
	QB_fprintf(outfile,"Position data Generated by SPaMD\n");
    QB_fprintf(outfile,"%d\tatoms\n",all);
    QB_fprintf(outfile,"%d\tatom types\n",QB->TypeNumber);
	double mat[3][3];
	double temp;
    if(QB->box)
	{
		//force transformation to lammps 6 par form
		double boundxx=sqrt(QB->mat[0][0]*QB->mat[0][0]+QB->mat[0][1]*QB->mat[0][1]+QB->mat[0][2]*QB->mat[0][2]);//found xx
		mat[0][0]=QB->mat[0][0]/boundxx;
		mat[0][1]=QB->mat[0][1]/boundxx;
		mat[0][2]=QB->mat[0][2]/boundxx;
		double boundxy=mat[0][0]*QB->mat[1][0]+mat[0][1]*QB->mat[1][1]+mat[0][2]*QB->mat[1][2];//vec.y*mat.x=vec.xy
		mat[1][0]=QB->mat[1][0]-boundxy*mat[0][0];
		mat[1][1]=QB->mat[1][1]-boundxy*mat[0][1];
		mat[1][2]=QB->mat[1][2]-boundxy*mat[0][2];
		double boundyy=sqrt(mat[1][0]*mat[1][0]+mat[1][1]*mat[1][1]+mat[1][2]*mat[1][2]);//|vec.y-vec.xy|
		mat[1][0]=mat[1][0]/boundyy;
		mat[1][1]=mat[1][1]/boundyy;
		mat[1][2]=mat[1][2]/boundyy;
		double boundxz=mat[0][0]*QB->mat[2][0]+mat[0][1]*QB->mat[2][1]+mat[0][2]*QB->mat[2][2];//vec.z*mat.x=vec.xz
		double boundyz=mat[1][0]*QB->mat[2][0]+mat[1][1]*QB->mat[2][1]+mat[1][2]*QB->mat[2][2];//vec.z*mat.y=vec.yz
		mat[2][0]=QB->mat[2][0]-boundxz*mat[0][0]-boundyz*mat[1][0];
		mat[2][1]=QB->mat[2][1]-boundxz*mat[0][1]-boundyz*mat[1][1];
		mat[2][2]=QB->mat[2][2]-boundxz*mat[0][2]-boundyz*mat[1][2];
		double boundzz=sqrt(mat[2][0]*mat[2][0]+mat[2][1]*mat[2][1]+mat[2][2]*mat[2][2]);//|vec.y-vec.xz-vec.yz|
		mat[2][0]=mat[2][0]/boundzz;
		mat[2][1]=mat[2][1]/boundzz;
		mat[2][2]=mat[2][2]/boundzz;
		QB_fprintf(outfile,"%lf\t\t%lf\txlo xhi\n",QB->zerox,QB->zerox+boundxx);
    	QB_fprintf(outfile,"%lf\t\t%lf\tylo yhi\n",QB->zeroy,QB->zeroy+boundyy);
    	QB_fprintf(outfile,"%lf\t\t%lf\tzlo zhi\n",QB->zeroz,QB->zeroz+boundzz);
		QB_fprintf(outfile,"%lf\t\t%lf\t%lf\txy xz yz\n",boundxy,boundxz,boundyz);
		QB_fprintf(outfile,"\nAtoms\n\n");
		double x,y,z;
    	for(i=0;i<QB->TotalNumber;i++)
    	{
    	    if(QB->atom[i].id>0)
    	    {
				x=(QB->atom[i].x-QB->zerox)*mat[0][0]+(QB->atom[i].y-QB->zeroy)*mat[0][1]+(QB->atom[i].z-QB->zeroz)*mat[0][2]+QB->zerox;
				y=(QB->atom[i].x-QB->zerox)*mat[1][0]+(QB->atom[i].y-QB->zeroy)*mat[1][1]+(QB->atom[i].z-QB->zeroz)*mat[1][2]+QB->zeroy;
				z=(QB->atom[i].x-QB->zerox)*mat[2][0]+(QB->atom[i].y-QB->zeroy)*mat[2][1]+(QB->atom[i].z-QB->zeroz)*mat[2][2]+QB->zeroz;
    	    	k++;
    	        QB_fprintf(outfile,"%d\t%d\t%lf\t%lf\t%lf\n",k,QB->atom[i].type,x,y,z);
    	    }
    	}
	}
	else
	{
		QB_fprintf(outfile,"%lf\t\t%lf\txlo xhi\n",QB->startx,QB->startx+QB->boundx);
    	QB_fprintf(outfile,"%lf\t\t%lf\tylo yhi\n",QB->starty,QB->starty+QB->boundy);
    	QB_fprintf(outfile,"%lf\t\t%lf\tzlo zhi\n",QB->startz,QB->startz+QB->boundz);
    	QB_fprintf(outfile,"\nAtoms\n\n");
    	for(i=0;i<QB->TotalNumber;i++)
    	{
    	    if(QB->atom[i].id>0)
    	    {
    	    	k++;
    	        QB_fprintf(outfile,"%d\t%d\t%lf\t%lf\t%lf\n",QB->atom[i].id,QB->atom[i].type,QB->atom[i].x,QB->atom[i].y,QB->atom[i].z);
    	    }
    	}
	}
    	QB_fclose(outfile);
}


void QB_dump_lmc(QB_tools *QB,const char *name)
{
    int i,all=0,j,k=0;
    FILE* outfile;
	outfile=QB_fopen(name,"w");
	for(i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].id>0)
			all++;
	}
	QB_fprintf(outfile,"ITEM: TIMESTEP\n0\n");
	QB_fprintf(outfile,"ITEM: NUMBER OF ATOMS\n%d\n",all);
	double mat[3][3];
	double temp;
	if(QB->box)
    {
		//force transformation to lammps 6 par form
		double boundxx=sqrt(QB->mat[0][0]*QB->mat[0][0]+QB->mat[0][1]*QB->mat[0][1]+QB->mat[0][2]*QB->mat[0][2]);//found xx
		mat[0][0]=QB->mat[0][0]/boundxx;
		mat[0][1]=QB->mat[0][1]/boundxx;
		mat[0][2]=QB->mat[0][2]/boundxx;
		double boundxy=mat[0][0]*QB->mat[1][0]+mat[0][1]*QB->mat[1][1]+mat[0][2]*QB->mat[1][2];//vec.y*mat.x=vec.xy
		mat[1][0]=QB->mat[1][0]-boundxy*mat[0][0];
		mat[1][1]=QB->mat[1][1]-boundxy*mat[0][1];
		mat[1][2]=QB->mat[1][2]-boundxy*mat[0][2];
		double boundyy=sqrt(mat[1][0]*mat[1][0]+mat[1][1]*mat[1][1]+mat[1][2]*mat[1][2]);//|vec.y-vec.xy|
		mat[1][0]=mat[1][0]/boundyy;
		mat[1][1]=mat[1][1]/boundyy;
		mat[1][2]=mat[1][2]/boundyy;
		double boundxz=mat[0][0]*QB->mat[2][0]+mat[0][1]*QB->mat[2][1]+mat[0][2]*QB->mat[2][2];//vec.z*mat.x=vec.xz
		double boundyz=mat[1][0]*QB->mat[2][0]+mat[1][1]*QB->mat[2][1]+mat[1][2]*QB->mat[2][2];//vec.z*mat.y=vec.yz
		mat[2][0]=QB->mat[2][0]-boundxz*mat[0][0]-boundyz*mat[1][0];
		mat[2][1]=QB->mat[2][1]-boundxz*mat[0][1]-boundyz*mat[1][1];
		mat[2][2]=QB->mat[2][2]-boundxz*mat[0][2]-boundyz*mat[1][2];
		double boundzz=sqrt(mat[2][0]*mat[2][0]+mat[2][1]*mat[2][1]+mat[2][2]*mat[2][2]);//|vec.y-vec.xz-vec.yz|
		mat[2][0]=mat[2][0]/boundzz;
		mat[2][1]=mat[2][1]/boundzz;
		mat[2][2]=mat[2][2]/boundzz;
		double startx,starty,startz;
		startx=QB->zerox;
		starty=QB->zeroy;
		startz=QB->zeroz;
		if(boundxy<0)startx+=boundxy;
		if(boundxz<0)startx+=boundxz;
		if(boundyz<0)starty+=boundyz;
		double boundx=boundxx+fabs(boundxy)+fabs(boundxz);
		double boundy=boundyy+fabs(boundyz);
		double boundz=boundzz;
		QB_fprintf(outfile,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n");	
		QB_fprintf(outfile,"%lf %lf %lf\n",startx,startx+boundx,boundxy);
		QB_fprintf(outfile,"%lf %lf %lf\n",starty,starty+boundy,boundxz);
		QB_fprintf(outfile,"%lf %lf %lf\n",startz,startz+boundz,boundyz);
		QB_fprintf(outfile,"ITEM: ATOMS id type x y z");
		for(i=0;i<QB->exchannel;i++)
        	QB_fprintf(outfile," %s",QB->exch[i]);
    	for(i=0;i<QB->num_exint;i++)
    		QB_fprintf(outfile," %s",QB->name_exint[i]);
        for(i=0;i<QB->num_exdouble;i++)
        	QB_fprintf(outfile," %s",QB->name_exdouble[i]);
        QB_fprintf(outfile,"\n");
		double x,y,z;
		for(i=0;i<QB->TotalNumber;i++)
    	{
    	    if(QB->atom[i].id>0)
    	    {
				x=(QB->atom[i].x-QB->zerox)*mat[0][0]+(QB->atom[i].y-QB->zeroy)*mat[0][1]+(QB->atom[i].z-QB->zeroz)*mat[0][2]+QB->zerox;
				y=(QB->atom[i].x-QB->zerox)*mat[1][0]+(QB->atom[i].y-QB->zeroy)*mat[1][1]+(QB->atom[i].z-QB->zeroz)*mat[1][2]+QB->zeroy;
				z=(QB->atom[i].x-QB->zerox)*mat[2][0]+(QB->atom[i].y-QB->zeroy)*mat[2][1]+(QB->atom[i].z-QB->zeroz)*mat[2][2]+QB->zeroz;
    	    	k++;
    	    	QB_fprintf(outfile,"%d\t%d\t%lf\t%lf\t%lf",k,QB->atom[i].type,x,y,z);
    	        for(j=0;j<QB->exchannel;j++)
    	        	QB_fprintf(outfile,"\t%lf",QB->exdat[j][i]);
    	        for(j=0;j<QB->num_exint;j++)
    	        	QB_fprintf(outfile,"\t%d",QB->par_exint[j][i]);
    	        for(j=0;j<QB->num_exdouble;j++)
    	        	QB_fprintf(outfile,"\t%lf",QB->par_exdouble[j][i]);
				QB_fprintf(outfile,"\n");
    	    }
    	}	
    }
    else
    {
        QB_fprintf(outfile,"ITEM: BOX BOUNDS pp pp pp\n");	
		QB_fprintf(outfile,"%lf %lf\n",QB->startx,QB->startx+QB->boundx);
		QB_fprintf(outfile,"%lf %lf\n",QB->starty,QB->starty+QB->boundy);
		QB_fprintf(outfile,"%lf %lf\n",QB->startz,QB->startz+QB->boundz);
		QB_fprintf(outfile,"ITEM: ATOMS id type x y z");
		for(i=0;i<QB->exchannel;i++)
        	QB_fprintf(outfile," %s",QB->exch[i]);
    	for(i=0;i<QB->num_exint;i++)
    		QB_fprintf(outfile," %s",QB->name_exint[i]);
        for(i=0;i<QB->num_exdouble;i++)
        	QB_fprintf(outfile," %s",QB->name_exdouble[i]);
        QB_fprintf(outfile,"\n");
    	for(i=0;i<QB->TotalNumber;i++)
    	{
    	    if(QB->atom[i].id>0)
    	    {
    	    	k++;
    	    	QB_fprintf(outfile,"%d\t%d\t%lf\t%lf\t%lf",QB->atom[i].id,QB->atom[i].type,QB->atom[i].x,QB->atom[i].y,QB->atom[i].z);
    	        for(j=0;j<QB->exchannel;j++)
    	        	QB_fprintf(outfile,"\t%lf",QB->exdat[j][i]);
    	        for(j=0;j<QB->num_exint;j++)
    	        	QB_fprintf(outfile,"\t%d",QB->par_exint[j][i]);
    	        for(j=0;j<QB->num_exdouble;j++)
    	        	QB_fprintf(outfile,"\t%lf",QB->par_exdouble[j][i]);
				QB_fprintf(outfile,"\n");
    	    }
    	}
	}
    QB_fclose(outfile);
}
#ifdef SPaMD_model
int QB_SPaMD_upload(QB_tools *QB)
{
    int i;
    vector r,v;
    v.x=v.y=v.z=0;
    int itag=0;
    SPaMD_memory(1.1*QB->TotalNumber/SPaMD_partition_size());
    SPaMD_geometry(QB->startx,QB->starty,QB->startz,QB->boundx,QB->boundy,QB->boundz,Cutoff);
    SPaMD_sync_with_nodes();
    for(i=0;i<QB->TotalNumber;i++)   
        if(QB->atom[i].id>0)
    	{	
   	    itag++;
    	    r.x=QB->atom[i].x;
    	    r.y=QB->atom[i].y;
    	    r.z=QB->atom[i].z;
	    save_particle(&r, &v, (QB->atom[i].type)-1, itag,0);
        }
    SPaMD_sync_with_nodes();
    int        np,i0,j0;
    Particle  *pt;
    i0=j0=0;
    pt=(Particle *) Cells[0][0][0].ptr;
    for(i=0;i<Maxsize;i++,pt++)
    {
        if (pt->type==0) i0++;
        if (pt->type==1) j0++;
    }

    i0=SPaMD_reduce_int(i0,SPaMD_combiner_add);
    j0=SPaMD_reduce_int(j0,SPaMD_combiner_add);

    SPaMD_fprintf(stdout,"Type 0: %d and Type 1: %d found\n", i0,j0);
    SPaMD_fprintf(stdout,"Total %d particles found\n", np = SPaMD_count_particles());

    /* Form the cells */
    SPaMD_form_cells();
    SPaMD_fprintf(stderr,"NEW BOX geometry created:\n");
    SPaMD_fprintf(stderr,"             %g %g\n", Pmin_x, Pmax_x);
    SPaMD_fprintf(stderr,"             %g %g\n", Pmin_y, Pmax_y);
    SPaMD_fprintf(stderr,"             %g %g\n", Pmin_z, Pmax_z);
    SPaMD_fprintf(stderr,"             %g %g\n", Bmin_x, Bmax_x);
    SPaMD_fprintf(stderr,"             %g %g\n", Bmin_y, Bmax_y);
    SPaMD_fprintf(stderr,"             %g %g\n", Bmin_z, Bmax_z);
	
    return np;
}
#endif

