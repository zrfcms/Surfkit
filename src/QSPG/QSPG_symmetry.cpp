#include"./../QB/QB.h"
#include"QSPG.h"

void QSPG_get_symmetry(QB_tools*input,QSPGDataset**data,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	if(spg.num!=0)
	{
		spg.pos=(double(*)[3])realloc(spg.pos,sizeof(double)*12*spg.num);
		spg.type=(int*)realloc(spg.type,sizeof(int)*4*spg.num);
	}
	*data=spg_get_dataset(spg.mat,
									 spg.pos,
									 spg.type,
									 spg.num,
									 symprec);
	if(spg.num!=0)
		SPG_free(&spg);
}

void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);

void SPG_add_symmetry_atom(SPG_tools *spg,int hall_number,int type,double x,double y,double z,double symprec);
void SPG_exchange_symmetry_atom(SPG_tools *spg,int hall_number,int type,double x,double y,double z,double symprec);
void QSPG_add_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,int type,double x,double y,double z,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	SPG_add_symmetry_atom(&spg,hall_number,type,x,y,z,symprec);
	SPG_wrap_box(&spg,symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
}

void QSPG_exchange_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,int type,double x,double y,double z,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	SPG_exchange_symmetry_atom(&spg,hall_number,type,x,y,z,symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
}

void SPG_add_symmetry_atom(SPG_tools *spg,int hall_number,int type,double x,double y,double z,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	int rot[192][3][3];
	double tra[192][3];
	int symmetry_num=spg_get_symmetry_from_database(rot,tra,hall_number);
	
	int failed_flag;
	double vec[3];
	double dsp[3];
  
	for(i=0;i<symmetry_num;i++)
	{
		failed_flag=0;
		vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
		vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
		vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
		for(j=0;j<spg->num;j++)
		{
			dsp[0]=vec[0]-spg->pos[j][0];
			dsp[1]=vec[1]-spg->pos[j][1];
			dsp[2]=vec[2]-spg->pos[j][2];
			while(dsp[0]<-0.5)dsp[0]+=1;
			while(dsp[1]<-0.5)dsp[1]+=1;
			while(dsp[2]<-0.5)dsp[2]+=1;
			while(dsp[0]>0.5)dsp[0]-=1;
			while(dsp[1]>0.5)dsp[1]-=1;
			while(dsp[2]>0.5)dsp[2]-=1;
			if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
			{
				failed_flag=1;
				spg->type[j]=type;
				break;
			}
		}
		if(!failed_flag)
			SPG_add_atom(spg,type,vec[0],vec[1],vec[2]);
	}
}

void SPG_exchange_symmetry_atom(SPG_tools *spg,int hall_number,int type,double x,double y,double z,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	int rot[192][3][3];
	double tra[192][3];
	int symmetry_num=spg_get_symmetry_from_database(rot,tra,hall_number);
	
	int failed_flag;
	double vec[3];
	double dsp[3];
  
	for(i=0;i<symmetry_num;i++)
	{
		vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
		vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
		vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
		for(j=0;j<spg->num;j++)
		{
			dsp[0]=vec[0]-spg->pos[j][0];
			dsp[1]=vec[1]-spg->pos[j][1];
			dsp[2]=vec[2]-spg->pos[j][2];
			while(dsp[0]<-0.5)dsp[0]+=1;
			while(dsp[1]<-0.5)dsp[1]+=1;
			while(dsp[2]<-0.5)dsp[2]+=1;
			while(dsp[0]>0.5)dsp[0]-=1;
			while(dsp[1]>0.5)dsp[1]-=1;
			while(dsp[2]>0.5)dsp[2]-=1;
			if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
			{
				spg->type[j]=type;
			}
		}
	}
}

void QSPG_add_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],int type,double x,double y,double z,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	
	int i,j;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	
	int failed_flag;
	double vec[3];
	double dsp[3];
  
	for(i=0;i<symmetry_num;i++)
	{
		failed_flag=0;
		vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
		vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
		vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
		for(j=0;j<spg.num;j++)
		{
			dsp[0]=vec[0]-spg.pos[j][0];
			dsp[1]=vec[1]-spg.pos[j][1];
			dsp[2]=vec[2]-spg.pos[j][2];
			while(dsp[0]<-0.5)dsp[0]+=1;
			while(dsp[1]<-0.5)dsp[1]+=1;
			while(dsp[2]<-0.5)dsp[2]+=1;
			while(dsp[0]>0.5)dsp[0]-=1;
			while(dsp[1]>0.5)dsp[1]-=1;
			while(dsp[2]>0.5)dsp[2]-=1;
			if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
			{
				failed_flag=1;
				spg.type[j]=type;
				break;
			}
		}
		if(!failed_flag)
			SPG_add_atom(&spg,type,vec[0],vec[1],vec[2]);
	}
	SPG_wrap_box(&spg,symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
}

void QSPG_exchange_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],int type,double x,double y,double z,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	
	int i,j;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	
	int failed_flag;
	double vec[3];
	double dsp[3];
  
	for(i=0;i<symmetry_num;i++)
	{
		failed_flag=0;
		vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
		vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
		vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
		for(j=0;j<spg.num;j++)
		{
			dsp[0]=vec[0]-spg.pos[j][0];
			dsp[1]=vec[1]-spg.pos[j][1];
			dsp[2]=vec[2]-spg.pos[j][2];
			while(dsp[0]<-0.5)dsp[0]+=1;
			while(dsp[1]<-0.5)dsp[1]+=1;
			while(dsp[2]<-0.5)dsp[2]+=1;
			while(dsp[0]>0.5)dsp[0]-=1;
			while(dsp[1]>0.5)dsp[1]-=1;
			while(dsp[2]>0.5)dsp[2]-=1;
			if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
			{
				spg.type[j]=type;
				break;
			}
		}
	}
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
}

void QSPG_simplify_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	
	int i,j,k;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	
	int failed_flag;
	double vec[3];
	double dsp[3];
	for(i=0;i<spg.num;i++)
	{
		if(spg.type[i]>0)
		for(j=i+1;j<spg.num;j++)
		{
			if(spg.type[j]>0)
			for(k=0;k<symmetry_num;k++)
			{
				vec[0]=rot[k][0][0]*spg.pos[i][0]+rot[k][0][1]*spg.pos[i][1]+rot[k][0][2]*spg.pos[i][2]+tra[k][0];
				vec[1]=rot[k][1][0]*spg.pos[i][0]+rot[k][1][1]*spg.pos[i][1]+rot[k][1][2]*spg.pos[i][2]+tra[k][1];
				vec[2]=rot[k][2][0]*spg.pos[i][0]+rot[k][2][1]*spg.pos[i][1]+rot[k][2][2]*spg.pos[i][2]+tra[k][2];
				dsp[0]=vec[0]-spg.pos[j][0];
				dsp[1]=vec[1]-spg.pos[j][1];
				dsp[2]=vec[2]-spg.pos[j][2];
				while(dsp[0]<-0.5)dsp[0]+=1;
				while(dsp[1]<-0.5)dsp[1]+=1;
				while(dsp[2]<-0.5)dsp[2]+=1;
				while(dsp[0]>0.5)dsp[0]-=1;
				while(dsp[1]>0.5)dsp[1]-=1;
				while(dsp[2]>0.5)dsp[2]-=1;
				if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
				{
					spg.type[j]=-1;
					break;
				}
			}
		}
	}
	for(i=0;i<spg.num;i++)
	{
		if(spg.type[i]>0)
		{
			double x,y,z;
			x=spg.pos[i][0];
			y=spg.pos[i][1];
			z=spg.pos[i][2];
			double temp=spg.pos[i][0]*spg.pos[i][0]+spg.pos[i][1]*spg.pos[i][1]+spg.pos[i][2]*spg.pos[i][2];
			for(k=0;k<symmetry_num;k++)
			{
				vec[0]=rot[k][0][0]*x+rot[k][0][1]*y+rot[k][0][2]*z+tra[k][0];
				vec[1]=rot[k][1][0]*x+rot[k][1][1]*y+rot[k][1][2]*z+tra[k][1];
				vec[2]=rot[k][2][0]*x+rot[k][2][1]*y+rot[k][2][2]*z+tra[k][2];
				//set the particle to the nearest to the center
				if(x>0)
				if(y>0)
				if(z>0)
				if(x*x+y*y+z*z<temp)
				{
					spg.pos[i][0]=vec[0];
					spg.pos[i][1]=vec[1];
					spg.pos[i][2]=vec[2];
				}
			}
		}
	}
	
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
	for(i=0;i<output->TotalNumber;i++)
	{
		
		if(output->atom[i].type<0)
			QB_delete_atom(output,i);
	}
	QB_affine_delete(output);
}

void QSPG_simplify_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,double symprec)
{
	int rot[192][3][3];
	double tra[192][3];
	int symmetry_num=spg_get_symmetry_from_database(rot,tra,hall_number);
	QSPG_simplify_operation_atom(input,output,symmetry_num,rot,tra,symprec);
}

void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z)
{
	if(spg->num==0)
	{
		spg->pos=(double(*)[3])malloc(sizeof(double)*3);
		spg->type=(int*)malloc(sizeof(int));
	}
	else
	{
		spg->pos=(double(*)[3])realloc(spg->pos,(spg->num+1)*sizeof(double)*3);
		spg->type=(int*)realloc(spg->type,(spg->num+1)*sizeof(int));
	}
	spg->pos[spg->num][0]=x;
	spg->pos[spg->num][1]=y;
	spg->pos[spg->num][2]=z;
	spg->type[spg->num]=type;
	spg->num+=1;
}