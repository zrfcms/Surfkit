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
void 	   QB_save_elename(QB_tools *QB,int i,char* name);
int QB_read_vasp(QB_tools *QB,const char name[1024])
{    
	QB_free_elename(QB);
    FILE    *fi2;
    double  x,y,z;
	int 	i,j;
    int     n,now;
    char    ch[256], cc[6]="Atoms";
	char    gc[256];
	char 	line[1024];
	char   *token;
	double  ratio;
    QB->px=QB->py=QB->pz=1;
	for(i=0;i<3;i++)for(j=0;j<3;j++)QB->mat[i][j]=0;
    fi2=QB_fopen(name,"r");
    QB_fscanf(fi2,"%*[^\n]%*c");
	
    if(QB_fscanf(fi2,"%lf ",&ratio)!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[0][0])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[0][1])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[0][2])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[1][0])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[1][1])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[1][2])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[2][0])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[2][1])!=1){QB_fclose(fi2);return -1;}
	if(QB_fscanf(fi2,"%lf ",&QB->mat[2][2])!=1){QB_fclose(fi2);return -1;}
	for(i=0;i<3;i++)for(j=0;j<3;j++)QB->mat[i][j]*=ratio;
	int typelist_n=0;
	int namelist_n=0;
	int *typelist=(int*)malloc(sizeof(int));
	int flag=0;
	
	QB_fgets(line,1024,fi2);
	token = strtok(line, " \t\n\r\f");
	while( token != NULL ) {
		QB_save_elename(QB,namelist_n,token);
		namelist_n++;
		token = strtok(NULL, " \t\n\r\f");
	}
	
	QB_fgets(line,1024,fi2);
	token = strtok(line, " \t\n\r\f");
	while( token != NULL ) {
		typelist=(int*)realloc(typelist,(typelist_n+1)*sizeof(int));
		typelist[typelist_n]=(int)QB_checkdat(token);
		typelist_n++;
		token = strtok(NULL, " \t\n\r\f");
	}
	QB_fgets(gc,1024,fi2);
	int fix_flag=0;
	if(gc[0]=='S'||gc[0]=='s')
	{
		fix_flag=1;
		QB_fgets(gc,1024,fi2);
	}
	while(gc[0]=='!'){QB_fscanf(fi2,"%*[^\n]%*c");if(QB_fscanf(fi2,"%s ",gc)!=1){QB_fclose(fi2);return -1;}}
	if(gc[0]=='C'||gc[0]=='c'||gc[0]=='K'||gc[0]=='k')flag=1;else flag=0;
	QB->TotalNumber=0;
	for(i=0;i<typelist_n;i++)
		QB->TotalNumber+=typelist[i];
    QB->atom=(atomdata*)malloc(QB->TotalNumber*sizeof(atomdata));
	
	int current_type=0;
	int fix_a_id=-1;
	int fix_b_id=-1;
	int fix_c_id=-1;
	if(fix_flag)
	{
		QB_add_exint(QB,"Fix_a");
		QB_add_exint(QB,"Fix_b");
		QB_add_exint(QB,"Fix_c");
		fix_a_id=QB_slot_get(QB,"Fix_a");
		fix_b_id=QB_slot_get(QB,"Fix_b");
		fix_c_id=QB_slot_get(QB,"Fix_c");
	}
	
	j=0;
	SV_wait("Reading model please wait");
	for(i=0;i<QB->TotalNumber;i++)
	{
		if((i%100000)==0)
			SV_wait_update(i*100/QB->TotalNumber);
		
		if(i==current_type){current_type+=typelist[j];j++;}
		
		QB_fgets(line,1024,fi2);
		token = strtok(line, " \t\n\r\f");
		if(token!=NULL)x=QB_checkdat(token);
		token = strtok(NULL, " \t\n\r\f");
		if(token!=NULL)y=QB_checkdat(token);
		token = strtok(NULL, " \t\n\r\f");
		if(token!=NULL)z=QB_checkdat(token);
		if(fix_flag)
		{
			token = strtok(NULL, " \t\n\r\f");
			if(token!=NULL)
			{
				if(token[0]=='F')QB_slot_save(QB,fix_a_id,i,1);
				else QB_slot_save(QB,fix_a_id,i,0);
			}
			token = strtok(NULL, " \t\n\r\f");
			if(token!=NULL)
			{
				if(token[0]=='F')QB_slot_save(QB,fix_b_id,i,1);
				else QB_slot_save(QB,fix_b_id,i,0);
			}
			token = strtok(NULL, " \t\n\r\f");
			if(token!=NULL)
			{
				if(token[0]=='F')QB_slot_save(QB,fix_c_id,i,1);
				else QB_slot_save(QB,fix_c_id,i,0);
			}
		}
		if(flag)
		{
			QB->atom[i].x=x*ratio;
			QB->atom[i].y=y*ratio;
			QB->atom[i].z=z*ratio;
		}
		else
		{
			QB->atom[i].x=x*QB->mat[0][0]+y*QB->mat[1][0]+z*QB->mat[2][0];
			QB->atom[i].y=x*QB->mat[0][1]+y*QB->mat[1][1]+z*QB->mat[2][1];
			QB->atom[i].z=x*QB->mat[0][2]+y*QB->mat[1][2]+z*QB->mat[2][2];
		}
		QB->atom[i].id=i+1;
		QB->atom[i].type=j;
	}
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
	QB->TypeNumber=typelist_n;
    QB_printf(">>VASP POSCAR file %s is found\n",name);
    QB_printf(">>total %d atoms\n",QB->TotalNumber);
    QB_printf(">>total %d atom types\n",QB->TypeNumber);
    QB_fclose(fi2);
	SV_wait_end();
    return 1;
}

void QB_dump_vasp_Cartesian(QB_tools *QB,const char *name)
{
	int i,j;
	FILE* outfile;
	double center[3];
	int fix_slot[3]={-1,-1,-1};
	int fix_flag=0;
	int fix_status[3];
	fix_slot[0]=QB_slot_get(QB,"Fix_a");
	fix_slot[1]=QB_slot_get(QB,"Fix_b");
	fix_slot[2]=QB_slot_get(QB,"Fix_c");
	outfile=QB_fopen(name,"w");
	QB_fprintf(outfile,"\tPosition data Generated by SPaMD\n");
    QB_fprintf(outfile,"1.0\n");
	if(QB->box)
	{
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n",QB->mat[0][0],QB->mat[0][1],QB->mat[0][2]);
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n",QB->mat[1][0],QB->mat[1][1],QB->mat[1][2]);
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n      ",QB->mat[2][0],QB->mat[2][1],QB->mat[2][2]);
		center[0]=QB->zerox;
		center[1]=QB->zeroy;
		center[2]=QB->zeroz;
	}
	else
	{
		QB_fprintf(outfile,"\t%lf\t\t0.0\t\t0.0\n",QB->boundx);
		QB_fprintf(outfile,"\t0.0\t\t%lf\t\t0.0\n",QB->boundy);
		QB_fprintf(outfile,"\t0.0\t\t0.0\t\t%lf\n      ",QB->boundz);
		center[0]=QB->startx;
		center[1]=QB->starty;
		center[2]=QB->startz;
	}
	char elename[256];
	int *typelist;
	typelist=(int*)malloc(QB->TypeNumber*sizeof(int));
	for(i=0;i<QB->TypeNumber;i++)
	{
		QB_load_elename(QB,i,elename);
		typelist[i]=0;		
	}
	for(i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].type>0&&QB->atom[i].type<=QB->TotalNumber)
		  if(QB->atom[i].id>0)
			typelist[QB_get_elename(QB,QB->ele.t[QB->atom[i].type-1].name)-1]++;
	}
	for(i=0;i<QB->TypeNumber;i++)
	if(typelist[i]>0)
	{
		QB_load_elename(QB,i,elename);
		QB_fprintf(outfile,"%s      ",elename);
	}
	QB_fprintf(outfile,"\n      ");
	for(i=0;i<QB->TypeNumber;i++)
		if(typelist[i]>0)
			QB_fprintf(outfile,"%d      ",typelist[i]);
	if((fix_slot[0]>0)&&(fix_slot[1]>0)&&(fix_slot[2]>0))
	{
		QB_fprintf(outfile,"\nSelective dynamics");
		fix_flag=1;
	}
	QB_fprintf(outfile,"\nCartesian\n");
	for(j=0;j<QB->TypeNumber;j++)
    for(i=0;i<QB->TotalNumber;i++)
    {
		if((QB_get_elename(QB,QB->ele.t[QB->atom[i].type-1].name)-1)==j)
        if(QB->atom[i].id>0)
		{
            QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf",QB->atom[i].x-center[0],QB->atom[i].y-center[1],QB->atom[i].z-center[2]);
			if(fix_flag)
			{
				fix_status[0]=(int)QB_slot_data(QB,fix_slot[0],i);
				fix_status[1]=(int)QB_slot_data(QB,fix_slot[1],i);
				fix_status[2]=(int)QB_slot_data(QB,fix_slot[2],i);
				if(fix_status[0]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
				if(fix_status[1]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
				if(fix_status[2]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
			}
			QB_fprintf(outfile,"\n");
		}
    }
	QB_fclose(outfile);
}

void QB_dump_vasp_Direct(QB_tools *QB,const char *name)
{
	int fix_slot[3]={-1,-1,-1};
	int fix_flag=0;
	int fix_status[3];
	fix_slot[0]=QB_slot_get(QB,"Fix_a");
	fix_slot[1]=QB_slot_get(QB,"Fix_b");
	fix_slot[2]=QB_slot_get(QB,"Fix_c");
	int i,j;
	FILE* outfile;
	double center[3];
	double temp[3];
	double d[3];
	double A;
	double mat[3][3];
	outfile=QB_fopen(name,"w");
	QB_fprintf(outfile,"\tPosition data Generated by SPaMD\n");
    QB_fprintf(outfile,"1.0\n");
	if(QB->box)
	{
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n",QB->mat[0][0],QB->mat[0][1],QB->mat[0][2]);
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n",QB->mat[1][0],QB->mat[1][1],QB->mat[1][2]);
		QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf\n      ",QB->mat[2][0],QB->mat[2][1],QB->mat[2][2]);
		center[0]=QB->zerox;
		center[1]=QB->zeroy;
		center[2]=QB->zeroz;
	}
	else
	{
		QB_fprintf(outfile,"\t%lf\t\t0.0\t\t0.0\n",QB->boundx);
		QB_fprintf(outfile,"\t0.0\t\t%lf\t\t0.0\n",QB->boundy);
		QB_fprintf(outfile,"\t0.0\t\t0.0\t\t%lf\n      ",QB->boundz);
		center[0]=QB->startx;
		center[1]=QB->starty;
		center[2]=QB->startz;
	}
	char elename[256];
	int *typelist;
	typelist=(int*)malloc(QB->TypeNumber*sizeof(int));
	for(i=0;i<QB->TypeNumber;i++)
	{
		QB_load_elename(QB,i,elename);//make sure every type get a name
		typelist[i]=0;
	}
	
	for(i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].type>0&&QB->atom[i].type<=QB->TotalNumber)
		  if(QB->atom[i].id>0)
			typelist[QB_get_elename(QB,QB->ele.t[QB->atom[i].type-1].name)-1]++;//in case of same element
	}
	//ignore type without atoms
	for(i=0;i<QB->TypeNumber;i++)
	{
		if(typelist[i]>0)
		{
			QB_load_elename(QB,i,elename);
			QB_fprintf(outfile,"%s      ",elename);
		}
	}
	QB_fprintf(outfile,"\n      ");
	
	for(i=0;i<QB->TypeNumber;i++)
		if(typelist[i]!=0)
			QB_fprintf(outfile,"%d      ",typelist[i]);
	if((fix_slot[0]>0)&&(fix_slot[1]>0)&&(fix_slot[2]>0))
	{
		QB_fprintf(outfile,"\nSelective dynamics");
		fix_flag=1;
	}	
	QB_fprintf(outfile,"\nDirect\n");
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
		temp[0]=1.0/(QB->boundx-QB->startx);
		temp[1]=1.0/(QB->boundy-QB->starty);
		temp[2]=1.0/(QB->boundz-QB->startz);
	}
	for(j=0;j<QB->TypeNumber;j++)
    for(i=0;i<QB->TotalNumber;i++)
    {
		
		if((QB_get_elename(QB,QB->ele.t[QB->atom[i].type-1].name)-1)==j)
        if(QB->atom[i].id>0)
		{
			if(QB->box)
			{
				d[0]=A*((QB->atom[i].x-center[0])*mat[0][0]+(QB->atom[i].y-center[1])*mat[0][1]+(QB->atom[i].z-center[2])*mat[0][2]);
				d[1]=A*((QB->atom[i].x-center[0])*mat[1][0]+(QB->atom[i].y-center[1])*mat[1][1]+(QB->atom[i].z-center[2])*mat[1][2]);
				d[2]=A*((QB->atom[i].x-center[0])*mat[2][0]+(QB->atom[i].y-center[1])*mat[2][1]+(QB->atom[i].z-center[2])*mat[2][2]);
			}
			else
			{
				d[0]=temp[0]*(QB->atom[i].x-QB->startx);
				d[1]=temp[1]*(QB->atom[i].y-QB->starty);
				d[2]=temp[2]*(QB->atom[i].z-QB->startz);
				
			}
			QB_fprintf(outfile,"\t%lf\t\t%lf\t\t%lf",d[0],d[1],d[2]);
			if(fix_flag)
			{
				fix_status[0]=QB_slot_data(QB,fix_slot[0],i);
				fix_status[1]=QB_slot_data(QB,fix_slot[1],i);
				fix_status[2]=QB_slot_data(QB,fix_slot[2],i);
				if(fix_status[0]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
				if(fix_status[1]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
				if(fix_status[2]==1)QB_fprintf(outfile,"\tF");else QB_fprintf(outfile,"\tT");
			}
			QB_fprintf(outfile,"\n");
		}
    }
	QB_fclose(outfile);
}

void QB_save_elename(QB_tools *QB,int i,char* name)
{
	char tempname[32];
	if(QB->ele.n==0)
	{
		QB->ele.t=(VASP_ELENAME*)malloc((i+1)*sizeof(VASP_ELENAME));
		for(int j=QB->ele.n;j<i;j++)
		{
			sprintf(QB->ele.t[j].name,"E%d",j);
		}
		QB->ele.n=i+1;
	}
	else if(QB->ele.n<i+1)
	{
		QB->ele.t=(VASP_ELENAME*)realloc(QB->ele.t,(i+1)*sizeof(VASP_ELENAME));
		for(int j=QB->ele.n;j<i;j++)
		{
			sprintf(QB->ele.t[j].name,"E%d",j);
		}
		QB->ele.n=i+1;
	}
	sprintf(QB->ele.t[i].name,"%s",name);
}

void QB_load_elename(QB_tools *QB,int i,char* name)
{
	char tempname[32];
	if(QB->ele.n==0)
	{
		QB->ele.t=(VASP_ELENAME*)malloc((i+1)*sizeof(VASP_ELENAME));
		for(int j=QB->ele.n;j<i+1;j++)
		{
			sprintf(QB->ele.t[j].name,"E%d",j);
		}
		QB->ele.n=i+1;
	}
	else if(QB->ele.n<i+1)
	{
		QB->ele.t=(VASP_ELENAME*)realloc(QB->ele.t,(i+1)*sizeof(VASP_ELENAME));
		for(int j=QB->ele.n;j<i+1;j++)
		{
			sprintf(QB->ele.t[j].name,"E%d",j);
		}
		QB->ele.n=i+1;
	}
	sprintf(name,"%s",QB->ele.t[i].name);
}

int QB_get_elename(QB_tools *QB,char*name)
{
	int i;
	
	if(QB->ele.n<QB->TypeNumber)
	{
		if(QB->ele.n==0)
			QB->ele.t=(VASP_ELENAME*)malloc(QB->TypeNumber*sizeof(VASP_ELENAME));
		else
			QB->ele.t=(VASP_ELENAME*)realloc(QB->ele.t,QB->TypeNumber*sizeof(VASP_ELENAME));
		for(int j=QB->ele.n;j<QB->TypeNumber;j++)
		{
			sprintf(QB->ele.t[j].name,"E%d",j);
		}
		QB->ele.n=QB->TypeNumber;
	}
	//return position if name exist
	for(i=0;i<QB->ele.n;i++)
	{
		if(strcmp(QB->ele.t[i].name,name)==0)
			return i+1;
	}
	//add name if name not exist
	if(QB->ele.n==0)
	{
		QB->ele.t=(VASP_ELENAME*)malloc(sizeof(VASP_ELENAME));	
		strcpy(QB->ele.t[0].name,name);
		QB->ele.n=1;
	}
	else
	{
		QB->ele.t=(VASP_ELENAME*)realloc(QB->ele.t,(QB->ele.n+1)*sizeof(VASP_ELENAME));
		strcpy(QB->ele.t[QB->ele.n].name,name);
		QB->ele.n+=1;
	}
	if(QB->ele.n>QB->TypeNumber)
		QB->TypeNumber=QB->ele.n;
	return QB->ele.n;
}
void QB_free_elename(QB_tools *QB)
{
	if(QB->ele.n!=0)
		free(QB->ele.t);
	QB->ele.n=0;
}
