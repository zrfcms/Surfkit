#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
#define MAX_MEMORY 1024
void SK_layerdata_save(char filename[1024],int list_n,SK_scaned_data* list)
{
	int i;
	char buffer[1024];
	FILE*f=fopen(filename,"w");
	fprintf(f,"Twist Angle,Rotate Angle1,Rotate Angle2,Sigma1,Sigma2,a Axial Mismatch,b Axial Mismatch2,Angle Mismatch,1M11,1M12,1M21,1M22,2M11,2M12,2M21,2M22\n");
	for(i=0;i<list_n;i++)
	{
		sprintf(buffer,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d,%d\n",
		list[i].theta1-list[i].theta2,
		list[i].theta1,list[i].theta2,list[i].sigma1,list[i].sigma2,
		list[i].mismatch[0],list[i].mismatch[1],list[i].mismatch[2],
		list[i].trans1[0][0],list[i].trans1[0][1],
		list[i].trans1[1][0],list[i].trans1[1][1],
		list[i].trans2[0][0],list[i].trans2[0][1],
		list[i].trans2[1][0],list[i].trans2[1][1]);
		fwrite(buffer,sizeof(char),strlen(buffer)*sizeof(char),f);
	}
	fclose(f);
}

void SK_layerdata_load(char filename[1024],int *list_n,SK_scaned_data** list)
{
	int i;
	char buffer[1024];
	FILE*f=fopen(filename,"r");
	double getd[7];
	int geti[8];
	char line[1024];
	char *token;
	int failed_flag;
	int N_list=(*list_n);
	fgets(line,1024,f);
	while(fgets(line,1024,f)!=NULL)
	{
		failed_flag=0;
		token=strtok(line, ",\t\n\r\f");
		for(int i=0;i<7;i++)
		{
			token=strtok(NULL, ",\t\n\r\f");
			if(token!=NULL)
				getd[i]=atof(token);
			else failed_flag=1;
		}
		
		for(int i=0;i<8;i++)
		{
			token=strtok(NULL, ",\t\n\r\f");
			if(token!=NULL)
				geti[i]=atoi(token);
			else failed_flag=1;
		}
		
		if(failed_flag)break;
		
		if(N_list==0)
			(*list)=(SK_scaned_data*)malloc(sizeof(SK_scaned_data));
		else 
			(*list)=(SK_scaned_data*)realloc((*list),(N_list+1)*sizeof(SK_scaned_data));
			
		(*list)[N_list].theta1=getd[0];
		(*list)[N_list].theta2=getd[1];
		(*list)[N_list].sigma1=getd[2];
		(*list)[N_list].sigma2=getd[3];
		
		(*list)[N_list].trans1[0][0]=geti[0];
		(*list)[N_list].trans1[0][1]=geti[1];
		(*list)[N_list].trans1[1][0]=geti[2];
		(*list)[N_list].trans1[1][1]=geti[3];
		(*list)[N_list].trans1[0][2]=(*list)[N_list].trans1[1][2]=
		(*list)[N_list].trans1[2][0]=(*list)[N_list].trans1[2][1]=0;
		(*list)[N_list].trans1[2][2]=1;
		
		(*list)[N_list].trans2[0][0]=geti[4];
		(*list)[N_list].trans2[0][1]=geti[5];
		(*list)[N_list].trans2[1][0]=geti[6];
		(*list)[N_list].trans2[1][1]=geti[7];
		(*list)[N_list].trans2[0][2]=(*list)[N_list].trans2[1][2]=
		(*list)[N_list].trans2[2][0]=(*list)[N_list].trans2[2][1]=0;
		(*list)[N_list].trans2[2][2]=1;
		
		(*list)[N_list].mismatch[0]=getd[4];
		(*list)[N_list].mismatch[1]=getd[5];
		(*list)[N_list].mismatch[2]=getd[6];
		(*list)[N_list].mismatch[3]=0;
		(*list)[N_list].mismatch[4]=0;
		
		(*list)[N_list].cylinder_flag=0;
		
		N_list++;
	}
	(*list_n)=N_list;
	fclose(f);
}