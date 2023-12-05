#include"./../QB/QB.h"
#include"QSPG.h"
int SPG_overloop_atom(SPG_tools* spg,int n,double symprec);
void QSPG_redefine(QB_tools* output,QB_tools* input,int rot[3][3],double symprec)
{
	int i,j;
	double mat[3][3];
	if((rot[0][0]*rot[1][1]*rot[2][2]+rot[0][1]*rot[1][2]*rot[2][0]+rot[0][2]*rot[1][0]*rot[2][1]-
	    rot[0][0]*rot[1][2]*rot[2][1]-rot[0][2]*rot[1][1]*rot[2][0]-rot[0][1]*rot[1][0]*rot[2][2])==0)
		return;
	if(input->TotalNumber==0)return;
	if(input!=output)//no atom added
	{
		QB_copy_system(output,input);
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mat[i][j]=input->mat[i][j];
	QB_duplicate(output,fabs(rot[0][0])+fabs(rot[1][0])+fabs(rot[2][0]),
					    fabs(rot[0][1])+fabs(rot[1][1])+fabs(rot[2][1]),
						fabs(rot[0][2])+fabs(rot[1][2])+fabs(rot[2][2]));
						
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			output->mat[i][j]=rot[i][0]*mat[0][j]+rot[i][1]*mat[1][j]+rot[i][2]*mat[2][j];					
	output->box=1;
	SPG_tools spg;
	spg.num=0;
	QB2SPG(output,&spg);
	for(i=spg.num-1;i>0;i--)
	{
		if(SPG_overloop_atom(&spg,i,symprec))spg.type[i]=-1;
	}
	SPG_wrap_box(&spg,symprec);
	SPG2QB(&spg,output);
	SPG_free(&spg);
	
	for(i=0;i<output->TotalNumber;i++)
	{
		if(output->atom[i].type<0)
			QB_delete_atom(output,i);
	}
	QB_affine_delete(output);
}

int SPG_overloop_atom(SPG_tools* spg,int n,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	double vec[3];	
	for(i=0;i<n;i++)
	{
		if(spg->type[i]!=spg->type[n])
			continue;
		vec[0]=spg->pos[i][0]-spg->pos[n][0];
		vec[1]=spg->pos[i][1]-spg->pos[n][1];
		vec[2]=spg->pos[i][2]-spg->pos[n][2];
		for(j=0;j<3;j++)
		{
			while(vec[j]<-0.5)vec[j]+=1;
			while(vec[j]>0.5)vec[j]-=1;
		}
		if((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])<symprec2)
		{
		   return 1;
		}
	}
	return 0;
}