#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
#define TIME_STEP 10
#define MAX_DISPLACE 0.1
double SK_rigid_average_r(QB_tools *QB,double r,double min,double max)
{
	QB->MCN=QB->TotalNumber;
	QB->pz=0; 
	QB_network_init(QB,r);
	double out=0;
	int num=0;
	double dis;
	for(int i=0;i<QB->TotalNumber;i++)
	{
		QB_list_build(QB,i,r);
		if(QB->atom[i].z>min&&QB->atom[i].z<=max)
		for(int j=0;j<QB->neb.num;j++)
			if(QB->atom[QB->neb.id[j]].z>min&&QB->atom[QB->neb.id[j]].z<=max)
			{

				dis=sqrt(QB_squaredlenths(QB->atom[i].x-QB->atom[QB->neb.id[j]].x,
										  QB->atom[i].y-QB->atom[QB->neb.id[j]].y,
										  QB->atom[i].z-QB->atom[QB->neb.id[j]].z));
				out+=dis;
				num++;
				break;
			}
		QB_list_clear(QB);
	}
	QB_network_free(QB);
	if(num)out/=num;
	return out;
}

static void SK_rigid_alias_move(QB_tools *QB,int slot,double d)
{
	QB->mat[2][2]+=d;
	for(int i=0;i<QB->TotalNumber;i++)
	{
		if(QB_slot_data(QB,slot,i))
			QB->atom[i].z+=TIME_STEP*d;
	}
	QB_bound_update(QB);
}

static void SK_rigid_alias_force_energy(QB_tools *QB,int slot,double max_cutoff,double r,double *force,double *energy)
{
	double d1=0,d2=0;
	int num=0;
	double dis;
	double rho=r/pow(2.0,1.0/6.0);
	for(int i=0;i<QB->MCN;i++)
	{
		if(!QB_slot_data(QB,slot,i))
		{
			QB_list_build(QB,i,max_cutoff);
			for(int j=0;j<QB->neb.num;j++)
			{
				if(QB_slot_data(QB,slot,QB->neb.id[j]))
				{
					dis=rho/sqrt(QB_squaredlenths(QB->atom[i].x-QB->atom[QB->neb.id[j]].x,
												  QB->atom[i].y-QB->atom[QB->neb.id[j]].y,
												  QB->atom[i].z-QB->atom[QB->neb.id[j]].z));
					d1+=(12*pow(dis,13)-6*pow(dis,7));
					d2+=(pow(dis,12)-pow(dis,6));
					num++;
					break;
				}
			}
			QB_list_clear(QB);
		}
	}
	if(num)
	{
		d1/=num;
		d2/=num;
	}
	(*force)=d1;
	(*energy)=d2;
}

static void SK_rigid_alias_energy(QB_tools *QB,int slot,double max_cutoff,double r,double *energy)
{
	double d1=0;
	int num=0;
	double dis;
	double rho=r/pow(2.0,1.0/6.0);
	for(int i=0;i<QB->MCN;i++)
	{
		if(!QB_slot_data(QB,slot,i))
		{
			QB_list_build(QB,i,max_cutoff);
			for(int j=0;j<QB->neb.num;j++)
			{
				if(QB_slot_data(QB,slot,QB->neb.id[j]))
				{
					dis=rho/sqrt(QB_squaredlenths(QB->atom[i].x-QB->atom[QB->neb.id[j]].x,
											      QB->atom[i].y-QB->atom[QB->neb.id[j]].y,
											      QB->atom[i].z-QB->atom[QB->neb.id[j]].z));
					d1+=(pow(dis,12)-pow(dis,6));
					num++;
					break;
				}
			}
			QB_list_clear(QB);
		}
	}
	if(num)
		d1/=num;
	(*energy)=d1;
}

void SK_rigid_alias(QB_tools *QB,double h,double r,int max_step,double symprec)
{
	QB->MCN=QB->TotalNumber;
	QB->pz=0;
	double max_cutoff=0;
	double max_z=QB->atom[0].z;
	double min_z=QB->atom[0].z;
	for(int i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].z>max_z)
			max_z=QB->atom[i].z;
		if(QB->atom[i].z<min_z)
			min_z=QB->atom[i].z;
	}
	max_cutoff=max_z-min_z;
	
	QB_pbc(QB,max_cutoff); 	
	QB->pz=1;
	QB_add_exint(QB,"Group");
	int slot=QB_slot_get(QB,"Group");
	QB_network_init(QB,2*r);
	
	for(int i=0;i<QB->TotalNumber;i++)
		QB_slot_save(QB,slot,i,(double)(QB->atom[i].z>h));
	
	int flag=0;
	int pair;
	double force;
	double energy1=0;
	double energy2=0;
	double rate=0.2;
	for(int i=0;i<max_step;i++)
	{
		SK_rigid_alias_force_energy(QB,slot,max_cutoff,r,&force,&energy1);
		printf("Step: %d;Force: %lf;Energy: %lf\n",i+1,force,energy1);
		if(fabs(force)>symprec)
		{
			for(int j=0;j<max_step;j++)
			{	
				double dis=force*pow(rate,j);
				if(fabs(dis)>MAX_DISPLACE)
					dis=fabs(dis)/dis*MAX_DISPLACE;
				SK_rigid_alias_move(QB,slot,dis);
				SK_rigid_alias_energy(QB,slot,max_cutoff,r,&energy2);
				if(energy2<energy1)break;
				if(fabs(force*pow(rate,j))<symprec){flag=1;break;}
				SK_rigid_alias_move(QB,slot,-dis);
			}
		}
		else break;
		if(flag)break;
	}
	QB_network_free(QB);
	QB_pbc_clean(QB);
}