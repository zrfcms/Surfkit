#include"./../QB/QB.h"
#include"QSPG.h"

void QSPG_primitive(QB_tools* output,QB_tools* input,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	spg.num=spg_find_primitive(spg.mat,
		      spg.pos,
		       spg.type,
		       spg.num,
		       symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no primitive cell found
	{
		QB_copy_system(output,input);
	}
}

void QSPG_refined(QB_tools* output,QB_tools* input,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	if(spg.num)
	{
		spg.pos=(double(*)[3])realloc(spg.pos,sizeof(double)*12*spg.num);
		spg.type=(int*)realloc(spg.type,4*sizeof(int)*spg.num);
	}
	spg.num=spg_refine_cell(spg.mat,
		      spg.pos,
		       spg.type,
		       spg.num,
		       symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no primitive cell found
	{
		QB_copy_system(output,input);
	}
}
//original methods get minimize cell of current scene
int SPG_overloop_vec(SPG_tools* spg,int type,double v_x,double v_y,double v_z,double symprec);
int QSPG_get_mintype(QB_tools* input);
extern int SPG_overloop_atom(SPG_tools* spg,int n,double symprec);
void QSPG_minimize(QB_tools* output,QB_tools* input,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	int min_type=QSPG_get_mintype(input);
	int start=0;
	int i,j,k,l;
	double mat[3]={1,1,1};
	
	for(i=0;i<spg.num;i++)
		if(spg.type[i]==min_type)
		{
			start=i;
			break;
		}
	
	double vec[3];
	double vec2[3];
	double temp;
	int failed_flag;
	for(i=0;i<spg.num;i++)
	{
		if(spg.type[i]==min_type)
		if(i!=start)
		{
			vec[0]=spg.pos[i][0]-spg.pos[start][0];
			vec[1]=spg.pos[i][1]-spg.pos[start][1];
			vec[2]=spg.pos[i][2]-spg.pos[start][2];
			//cos theta == 1 == parallel
			for(j=0;j<3;j++)
			{
				if(fabs(vec[0])+fabs(vec[1])+fabs(vec[2])-fabs(vec[j])<2*symprec)
				{
					
					failed_flag=0;
					for(k=0;k<spg.num;k++)
					{
						vec2[0]=spg.pos[k][0]+vec[0];
						vec2[1]=spg.pos[k][1]+vec[1];
						vec2[2]=spg.pos[k][2]+vec[2];
				
						if(!SPG_overloop_vec(&spg,spg.type[k],vec2[0],vec2[1],vec2[2],symprec))
						{
							failed_flag=1;
							break;
						}
					}
					if(!failed_flag)
						if(fabs(mat[j])>fabs(vec[j]))
							mat[j]=fabs(vec[j]);
				}
			}
		}
	}
	
	for(i=0;i<spg.num;i++)
	{
		for(j=0;j<3;j++)
			spg.pos[i][j]/=mat[j];
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			spg.mat[i][j]*=mat[j];
	
	for(i=spg.num-1;i>0;i--)
	{
		if(SPG_overloop_atom(&spg,i,symprec))spg.type[i]=-1;
	}
	SPG_wrap_box(&spg,symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no primitive cell found
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

int SPG_overloop_vec(SPG_tools* spg,int type,double v_x,double v_y,double v_z,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	double vec[3];	
	for(i=0;i<spg->num;i++)
	{
		if(spg->type[i]!=type)
			continue;
		vec[0]=spg->pos[i][0]-v_x;
		vec[1]=spg->pos[i][1]-v_y;
		vec[2]=spg->pos[i][2]-v_z;
		for(j=0;j<3;j++)
		{
			if(vec[j]<-0.5)vec[j]+=1;
			if(vec[j]>0.5)vec[j]-=1;
		}
		if((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])<symprec2)
		{
		   return 1;
		}
	}
	return 0;
}

int QSPG_get_mintype(QB_tools* input)
{
	int i;
	int ct=input->TotalNumber;
	int rt=1;
	//for safety check TypeNum first
	for(i=0;i<input->TotalNumber;i++)
	{
		if(input->TypeNumber<input->atom[i].type)
			input->TypeNumber=input->atom[i].type;
	}
	if(input->TypeNumber==1)return 1;
	int *type_list=(int*)malloc(sizeof(int)*input->TypeNumber);
	for(i=0;i<input->TypeNumber;i++)
		type_list[i]=0;
	
	for(i=0;i<input->TotalNumber;i++)
	{
		if(input->atom[i].type>0)
			type_list[input->atom[i].type-1]++;
	}
	for(i=0;i<input->TypeNumber;i++)
		if(type_list[i]!=0)
			if(type_list[i]<ct)
			{
				ct=type_list[i];
				rt=i+1;
			}
	free(type_list);
    return rt;
}

void SPG_wrap_box(SPG_tools* spg,double symprec)
{
	int i,j;
	for(i=0;i<spg->num;i++)
	for(j=0;j<3;j++)
	{
		if (spg->pos[i][j] < -symprec)
			spg->pos[i][j] = spg->pos[i][j] + 1.0 - (int) spg->pos[i][j];
		else
			spg->pos[i][j] = spg->pos[i][j] - (int) spg->pos[i][j];
	}
}